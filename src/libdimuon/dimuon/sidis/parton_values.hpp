#ifndef DM_SIDIS_PARTON_VALUES_H
#define DM_SIDIS_PARTON_VALUES_H

#include <cstddef>
#include <ranges>
#include <utility>

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/decay/decay_function.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/sidis/decay_parameters.hpp"
#include "dimuon/sidis/coefficient_channel.hpp"
#include "dimuon/sidis/distribution_matrix.hpp"
#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/sidis/parameters.hpp"
#include "dimuon/sidis/vectorization.hpp"

#include "dimuon/utility/math.hpp"
#include "dimuon/utility/shims.hpp"

namespace dm {

namespace sidis {
	/// @brief PDF and FF values at one integration point, shared by all partonic channels and structure functions
	/// @details Multiple channels need the same values, so each distinct value is computed once per integration point.
	/// The values are computed for all flavors at once (note that interpolation kernels are flavor-independent).
	/// The FF values are summed over the FF sets, each weighted by the decay function of its decay channel.
	/// The coefficient functions are linear in the FF values, so a single evaluation with these sums gives the sum over the decay channels.
	/// @tparam pdf_interface_t Type of PDF object. Must satisfy the `is_pdf_interface` concept.
	/// @tparam ff_interface_t Type of FF object. Must satisfy the `is_pdf_interface` concept.
	/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`). Must satisfy the `is_coupling` concept.
	template <is_pdf_interface pdf_interface_t, is_pdf_interface ff_interface_t, is_coupling coupling_t>
	class parton_values {
	public:
		/// @brief Computes the PDF and FF values at the given integration point
		/// @param[in] xi The PDF-side convolution variable
		/// @param[in] zeta The FF-side convolution variable
		/// @param[in] z The FF-side momentum fraction
		/// @param[in] params_ Parameters used for the evaluation of the SIDIS structure functions, with the setup of the integral, which must outlive the values
		/// @param[in] decay_params Parameters of the decay channels, one for each FF set, which give the weights of the FF sets
		template <is_decay_function decay_function_t>
		explicit parton_values(
			const double xi, const double zeta, const double z,
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params_,
			const decay_parameters<decay_function_t> &decay_params
		) :
			params{params_},
			pdf_matrices{compute_pdfs(xi, params_)},
			ff_matrices{compute_ffs(xi, zeta, z, params_, decay_params)} {}

		/// @brief Returns the precomputed PDF values of \p channel as a sparse vector of the rows of a matrix (see `sidis::distribution_matrix`). 
		/// The output contains all available PDF values at `momentum_fraction` and `momentum_fraction / xi` (`pdf` and `pdf_hat`).
		/// @param[in] channel The partonic channel (must be one of the channels used in the constructor)
		/// @param[in] coefficient_channel The coefficient channel of \p channel, which gives the PDF momentum fraction
		decltype(auto) pdf_rows(const sidis::partonic_channel &channel, const sidis::coefficient_channel &coefficient_channel) const {
			const std::size_t index = params.momentum_fraction_index(coefficient_channel.pdf_momentum_fraction(params.kinematics.x, params.kinematics.Q2));
			return pdf_matrices.at(index).rows(channel.initial_line.flavor);
		}

		/// @brief Returns the precomputed FF values of \p channel as a sparse vector of the rows of a matrix (see `sidis::distribution_matrix`).
		/// The output contains all available FF values at `z`, `z / zeta`, and `z / zeta_p` (`ff`, `ff_hat`, and `ff_hat_p`).
		/// @param[in] channel The partonic channel (must be one of the channels used in the constructor)
		/// @param[in] coefficient_channel The coefficient channel of \p channel, which gives the mass of the fragmenting parton
		decltype(auto) ff_rows(const sidis::partonic_channel &channel, const sidis::coefficient_channel &coefficient_channel) const {
			const std::size_t index = params.fragmentation_mass_index(coefficient_channel.fragmentation_mass);
			return ff_matrices.at(index).rows(channel.fragmentation_line.flavor);
		}

		/// @brief Computes the FF-side convolution variable `zeta_p`, which is related to `zeta` through the mass of the fragmenting parton.
		/// For a massless fragmenting parton, `zeta_p == zeta`.
		/// @param[in] xi The PDF-side convolution variable
		/// @param[in] zeta The FF-side convolution variable
		/// @param[in] Q2 The virtuality
		/// @param[in] fragmentation_mass The mass of the fragmenting parton
		static double zeta_p(const double xi, const double zeta, const double Q2, const double fragmentation_mass) noexcept {
			const double lambda = Q2 / (Q2 + math::power<2>(fragmentation_mass));
			return lambda == 1.0 ? zeta : ((1.0 - lambda) * xi + (1.0 - xi) * zeta) / (1.0 - lambda * xi);
		}

	private:
		/// The matrix of the PDF values at `momentum_fraction` and `momentum_fraction / xi` (`pdf` and `pdf_hat`)
		using pdf_matrix = sidis::distribution_matrix::type<pdf_interface_t, sidis::vectorization::pdf_value_count>;
		/// The matrix of the FF values at `z`, `z / zeta`, and `z / zeta_p` (`ff`, `ff_hat`, and `ff_hat_p`)
		using ff_matrix = sidis::distribution_matrix::type<ff_interface_t, sidis::vectorization::ff_value_count>;
		/// The largest number of distinct PDF momentum fractions of the channels, for which the values have room
		static constexpr std::size_t max_momentum_fractions = 36;
		/// The largest number of distinct fragmenting-parton masses of the channels, for which the values have room
		static constexpr std::size_t max_fragmentation_masses = 6;
		/// The PDF matrices, one for each PDF momentum fraction, in place (more throw `std::bad_alloc`)
		using pdf_matrices_storage = shims::inplace_vector<pdf_matrix, max_momentum_fractions>;
		/// The FF matrices, one for each fragmenting-parton mass, in place (more throw `std::bad_alloc`)
		using ff_matrices_storage = shims::inplace_vector<ff_matrix, max_fragmentation_masses>;

		/// The parameters of the integral
		const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params;
		/// The PDF matrices, in the order of `params.momentum_fractions`
		const pdf_matrices_storage pdf_matrices;
		/// The FF matrices, in the order of `params.fragmentation_masses`
		const ff_matrices_storage ff_matrices;

		/// Returns the PDF values at `momentum_fraction` (cached), and at `momentum_fraction / xi` (calculated here)
		static pdf_matrices_storage compute_pdfs(
			const double xi, 
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params_
		) {
			pdf_matrices_storage matrices;
			for (const auto &[momentum_fraction, pdfs_at_momentum_fractions] : std::views::zip(params_.momentum_fractions, params_.pdf_values_at_momentum_fractions)) {
				matrices.emplace_back(
					// PDF values at the momentum fraction, which is independent of the integration point and has therefore been precomputed
					pdfs_at_momentum_fractions,
					// PDF values at (momentum fraction) / xi, which depends on the integration point and must be calculated for every point
					distribution_matrix::compute_values(params_.pdf, momentum_fraction / xi, params_.factorization_scale)
				);
			}
			return matrices;
		}

		/// Returns the FF values at `z`, `z / zeta`, and `z / zeta_p` for each fragmenting-parton mass, summed over the FF sets.
		/// Each contribution to the FF sum is weighted by the decay function of its decay channel.
		template <is_decay_function decay_function_t>
		static ff_matrices_storage compute_ffs(
			const double xi, const double zeta, const double z, 
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params_, 
			const decay_parameters<decay_function_t> &decay_params
		) {
			// For the arithmetic operators on std::array
			using collections::operator+=, collections::operator*;
			
			// FF output type
			using ff_flavor_values = decltype(std::declval<const ff_interface_t &>().xfxQ2(0.0, 0.0));

			// Sum of FF values at `z` and `z / zeta`, weighted by the decay function
			ff_flavor_values ff_total{};
			ff_flavor_values ff_hat_total{};
			// Sum of FF values at `z / zeta_p`, weighted by the decay function and stored per the mass of the fragmenting parton (`zeta_p` depends on the mass)
			shims::inplace_vector<ff_flavor_values, max_fragmentation_masses> ff_hat_p_totals(params_.fragmentation_masses.size());

			for (const auto &[ff, decay] : std::views::zip(params_.ffs, decay_params.decay_instances)) {
				// Cut on the z-integral, which depends on the decay channel
				if (z < decay.z_min(params_.kinematics)) { continue; }

				// Skip decay channels that don't contribute
				const double weight = decay(params_.kinematics.x, z, params_.kinematics.Q2);
				if (weight == 0.0) { continue; }

				// Contributions to the totals from this decay channel
				ff_total += weight * distribution_matrix::compute_values(ff, z, params_.fragmentation_scale);
				ff_hat_total += weight * distribution_matrix::compute_values(ff, z / zeta, params_.fragmentation_scale);

				for (auto &&[ff_hat_p_total, fragmentation_mass] : std::views::zip(ff_hat_p_totals, params_.fragmentation_masses)) {
					ff_hat_p_total += weight * distribution_matrix::compute_values(ff, z / zeta_p(xi, zeta, params_.kinematics.Q2, fragmentation_mass), params_.fragmentation_scale);
				}
			}

			ff_matrices_storage matrices;
			for (const ff_flavor_values &ff_hat_p_total : ff_hat_p_totals) {
				matrices.emplace_back(ff_total, ff_hat_total, ff_hat_p_total);
			}
			return matrices;
		}
	};
}

}

#endif
