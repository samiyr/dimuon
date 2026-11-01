#ifndef DM_SIDIS_PARAMETERS_H
#define DM_SIDIS_PARAMETERS_H

#include <vector>
#include <flat_set>

#include "dimuon/common/coefficient_implementation.hpp"
#include "dimuon/common/heavy_quark_thresholds.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/coupling/is_coupling.hpp"
#include "dimuon/pdf/is_pdf_interface.hpp"
#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/sidis/coefficient_channel.hpp"
#include "dimuon/sidis/distribution_matrix.hpp"

namespace dm {

namespace sidis {
	/// Container for holding parameters related to the entire SIDIS cross section calculation
	template <is_pdf_interface pdf_interface_t, is_pdf_interface ff_interface_t, is_coupling coupling_t>
	struct parameters {
		/// PDF
		pdf_interface_t pdf;
		/// List of FFs
		std::vector<ff_interface_t> ffs;
		/// Coupling
		coupling_t coupling;
		/// Scattering process
		scattering_process process;
		/// Kinematical point
		trf_kinematics kinematics;
		/// The perturbative order of the channels and coefficient functions
		perturbative_order pto;
		/// The implementation of the coefficient functions at `pto`
		coefficient_implementation_t coefficient_implementation;
		/// Heavy-quark thresholds, which determine the number of active flavors
		heavy_quark_thresholds thresholds;
		/// Renormalization scale (GeV^2) at the kinematical point `kinematics`
		double renormalization_scale;
		/// Factorization scale (GeV^2) at the kinematical point `kinematics`
		double factorization_scale;
		/// Fragmentation scale (GeV^2) at the kinematical point `kinematics`
		double fragmentation_scale;
		/// Set of unique PDF-side momentum fractions (i.e. x and all possible values of the slow-rescaling variable across all channels)
		std::flat_set<double> momentum_fractions;
		/// A list of PDF values at `momentum_fractions`. The value type of the vector is the same as the output type of the PDF.
		std::vector<decltype(std::declval<const pdf_interface_t &>().xfxQ2(0.0, 0.0))> pdf_values_at_momentum_fractions;
		/// Set of unique fragmentation masses across all channels
		std::flat_set<double> fragmentation_masses;

		explicit parameters(
			pdf_interface_t pdf_,
			std::vector<ff_interface_t> ffs_,
			coupling_t coupling_,
			scattering_process process_,
			trf_kinematics kinematics_,
			perturbative_order pto_,
			const std::vector<sidis::partonic_channel> &channels,
			coefficient_implementation_t coefficient_implementation_,
			heavy_quark_thresholds thresholds_,
			double renormalization_scale_,
			double factorization_scale_,
			double fragmentation_scale_
		) :
		pdf{std::move(pdf_)},
		ffs{std::move(ffs_)},
		coupling{std::move(coupling_)},
		process{std::move(process_)},
		kinematics{std::move(kinematics_)},
		pto{pto_},
		coefficient_implementation{coefficient_implementation_},
		thresholds{std::move(thresholds_)},
		renormalization_scale{renormalization_scale_},
		factorization_scale{factorization_scale_},
		fragmentation_scale{fragmentation_scale_},
		momentum_fractions{
			channels 
			| std::views::transform([&](const sidis::partonic_channel &channel) {
				return sidis::coefficient_channel{channel}.pdf_momentum_fraction(kinematics.x, kinematics.Q2);
			})
			| std::ranges::to<std::flat_set<double>>()
		},
		pdf_values_at_momentum_fractions{
			momentum_fractions
			| std::views::transform([&](const double momentum_fraction) { 
				return distribution_matrix::compute_values(pdf, momentum_fraction, factorization_scale); 
			})
			| std::ranges::to<std::vector>()
		},
		fragmentation_masses{
			channels
			| std::views::transform([](const sidis::partonic_channel &channel) { return channel.fragmentation_line.mass; })
			| std::ranges::to<std::flat_set<double>>()
		} {}

		/// @brief Returns the position of \p momentum_fraction in `momentum_fractions`
		/// @details It is assumed that \p momentum_fraction is contained in `momentum_fractions` 
		std::size_t momentum_fraction_index(const double momentum_fraction) const noexcept {
			return collections::iterator_to_index(momentum_fractions, momentum_fractions.find(momentum_fraction));
		}

		/// @brief Returns the position of \p fragmentation_mass in `fragmentation_masses`
		/// @details It is assumed that \p fragmentation_mass is contained in `fragmentation_masses`
		std::size_t fragmentation_mass_index(const double fragmentation_mass) const noexcept {
			return collections::iterator_to_index(fragmentation_masses, fragmentation_masses.find(fragmentation_mass));
		}

	};
}

}

#endif