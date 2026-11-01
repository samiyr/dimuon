#ifndef DM_SIDIS_EVALUATOR_H
#define DM_SIDIS_EVALUATOR_H

#include <array>
#include <cassert>
#include <cmath>
#include <flat_map>
#include <flat_set>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#include "dimuon/common/coefficient_implementation.hpp"
#include "dimuon/common/is_observable.hpp"
#include "dimuon/common/partonic_channel_group.hpp"
#include "dimuon/common/perturbative_order.hpp"

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/decay/decay_function.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/sidis/coefficients/coefficient_function.hpp"
#include "dimuon/sidis/coefficients/evaluation_cache_gc.hpp"

#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/sidis/coefficient_channel.hpp"
#include "dimuon/sidis/physical_channel.hpp"
#include "dimuon/sidis/decay_parameters.hpp"
#include "dimuon/sidis/parameters.hpp"
#include "dimuon/sidis/parton_values.hpp"
#include "dimuon/sidis/vectorization.hpp"

#include "dimuon/utility/collections.hpp"
#include "dimuon/utility/vectorization.hpp"


namespace dm {

using distribution_index = std::size_t;

namespace sidis {
	/// A component filler is called from `sidis::evaluator::evaluate_observable` with
	/// - the PDF row index
	/// - the FF row index
	/// - the physical channel index
	/// - the actual scalar contribution
	/// and is responsible for filling the appropriate integral component with the value
	template <typename function_t>
	concept component_filler = requires(
		function_t fill_component, const distribution_index pdf_row, const distribution_index ff_row, const channel_index physical_channel_index, const double value
	) {
		fill_component(pdf_row, ff_row, physical_channel_index, value);
	};

	/// Utility class for computing SIDIS cross section values. 
	/// This class is responsible for actually computing integrands.
	/// For more details, see the `numerical_evaluation.md` document.
	class evaluator {
	public:
		/// @brief Constructs the evaluator, groups the partonic channels by their coefficient functions, and indexes each by the position of its physical channel in \p physical_channels
		/// @param[in] channels The partonic channels
		/// @param[in] physical_channels The physical channels (see `sidis::physical_channel`), which must include those of \p channels
		explicit evaluator(const std::vector<sidis::partonic_channel> &channels, const std::flat_set<sidis::physical_channel> &physical_channels) :
			coefficient_channels(group_by_coefficients(channels, physical_channels)) {}

		/// @brief Evaluates the coefficients of the given coefficient function (see `sidis::vectorization::distribution_product`)
		/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_pdf_interface` concept.
		/// @tparam ff_interface_t Type of FF object, used for obtaining FF values. Must satisfy the `is_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @param[in] x The PDF-side momentum fraction
		/// @param[in] z The FF-side momentum fraction
		/// @param[in] xi The PDF-side convolution variable
		/// @param[in] zeta The FF-side convolution variable
		/// @param[in] params Parameters used for the evaluation of the SIDIS structure functions
		/// @param[in] integrand The coefficient function of the partonic channels
		/// @param[in] channel The coefficient channel of the partonic channels (see `sidis::coefficient_channel`)
		/// @param[in] evaluation_cache Cached per-integration-point values that do not depend on the channel and are shared among the coefficient functions
		/// @param[in] Q2 The virtuality (GeV^2)
		/// @return The coefficient of each product of a PDF value and an FF value in the coefficient function (does not include couplings)
		template <is_pdf_interface pdf_interface_t, is_pdf_interface ff_interface_t, is_coupling coupling_t>
		sidis::vectorization::distribution_product evaluate_integrand(
			const double x, 
			const double z, 
			const double xi, 
			const double zeta, 
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params, 
			sidis::coefficients::coefficient_function auto &&integrand, 
			const sidis::coefficient_channel &channel,
			const sidis::evaluation_cache_gc &evaluation_cache,
			const double Q2
		) const {
			// Decide the PDF momentum fraction variable based on the mass scheme of the channels
			const double pdf_momentum_fraction = channel.pdf_momentum_fraction(x, Q2);

			// Enforce the lower bound of the xi-integral.
			// No check for zeta < ff_momentum_fraction, because that's taken into account elsewhere.
			// For the hw coefficient functions, the check is either inside the coefficient functions
			// or the coefficient functions are actually written as integrals from 0 to 1 and contain explicit theta function.
			// For the gc coefficient functions, the check is in `evaluate_observable`.
			if (xi < pdf_momentum_fraction) { return 0.0; }

			// Compute the coefficients by evaluating the integrand with unit PDF and FF values (see `sidis::vectorization::make_masking_vector`)
			using enum sidis::vectorization::distribution_value;
			return std::forward<decltype(integrand)>(integrand)({
				.channel = channel,
				.c = evaluation_cache,
				.xi = xi, 
				.pdf_momentum_fraction = pdf_momentum_fraction, 
				.pdf = sidis::vectorization::make_masking_vector(pdf),
				.pdf_hat = sidis::vectorization::make_masking_vector(pdf_hat),
				.zeta = zeta,
				.z = z,
				.ff = sidis::vectorization::make_masking_vector(ff),
				.ff_hat = sidis::vectorization::make_masking_vector(ff_hat),
				.ff_hat_p = sidis::vectorization::make_masking_vector(ff_hat_p),
				.factorization_scale = params.factorization_scale, 
				.fragmentation_scale = params.fragmentation_scale,
				.Q2 = Q2
			});
		}

		/// @brief Evaluates the coefficients of the given structure function for the partonic channels with a given coefficient channel (see `sidis::vectorization::distribution_product`)
		/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_pdf_interface` concept.
		/// @tparam ff_interface_t Type of FF object, used for obtaining FF values. Must satisfy the `is_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @param[in] xi The PDF-side convolution variable
		/// @param[in] zeta The FF-side convolution variable
		/// @param[in] z The FF-side momentum fraction
		/// @param[in] params Parameters used for the evaluation of the SIDIS structure functions
		/// @param[in] quark_to_quark The quark-to-quark structure function
		/// @param[in] quark_to_gluon The quark-to-gluon structure function
		/// @param[in] gluon_to_quark The gluon-to-quark structure function
		/// @param[in] gluon_to_gluon The gluon-to-gluon structure function
		/// @param[in] channel The coefficient channel of the partonic channels (see `sidis::coefficient_channel`)
		/// @param[in] evaluation_cache Values of the coefficient functions that do not depend on the channel, at this integration point (see `sidis::evaluation_cache_gc`)
		/// @return The coefficient of each product of a PDF value and an FF value (does not include couplings), which is the same for all the channels
		template <is_pdf_interface pdf_interface_t, is_pdf_interface ff_interface_t, is_coupling coupling_t>
		sidis::vectorization::distribution_product evaluate_coefficients(
			const double xi, const double zeta, const double z, 
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &params, 
			sidis::coefficients::coefficient_function auto &&quark_to_quark, 
			sidis::coefficients::coefficient_function auto &&quark_to_gluon,
			sidis::coefficients::coefficient_function auto &&gluon_to_quark, 
			sidis::coefficients::coefficient_function auto &&gluon_to_gluon,
			const sidis::coefficient_channel &channel,
			const sidis::evaluation_cache_gc &evaluation_cache
		) const {
			// Bound function for evaluating the given integrand
			const auto compute = [&](auto &&integrand) {
				return evaluate_integrand(
					params.kinematics.x, z, xi, zeta,
					params, 
					std::forward<decltype(integrand)>(integrand), 
					channel, evaluation_cache, params.kinematics.Q2
				);
			};

			// Compute the coefficients, deciding on which provided integrand function to
			// use based on the channel type.
			// For NNLO corrections, further divisions between different channel types will be
			// determined in the coefficient function itself.
			if (channel.is_quark_to_quark()) {
				return compute(std::forward<decltype(quark_to_quark)>(quark_to_quark));
			} else if (channel.is_quark_to_gluon()) {
				return compute(std::forward<decltype(quark_to_gluon)>(quark_to_gluon));
			} else if (channel.is_gluon_to_quark()) {
				return compute(std::forward<decltype(gluon_to_quark)>(gluon_to_quark));
			} else if (channel.is_gluon_to_gluon()) {
				return compute(std::forward<decltype(gluon_to_gluon)>(gluon_to_gluon));
			}
			std::unreachable();
		}

		/// @brief Evaluates the given observable at an integration point, as a sum of contributions: one for each partonic channel and each pair of a row of its PDF values
		/// and a row of its FF values
		///
		/// The PDF values of a channel are a matrix with a row for each function that takes the place of the PDF and a column for each PDF value of the coefficient functions
		/// (`pdf` and `pdf_hat`): a single row for a PDF whose values are numbers, such as a PDF set, or a row for each function of one whose values are vectors,
		/// such as the interpolation kernels (see `is_pdf_interface` and `sidis::distribution_matrix`).
		/// The FF values are a matrix in the same way, with the columns `ff`, `ff_hat`, and `ff_hat_p`. The matrices are sparse vectors of their rows (see `sparse_vector`):
		/// only the non-zero rows are stored.
		/// The coefficients of each group of channels are a matrix with a row for each PDF value and a column for each FF value (see `sidis::vectorization::distribution_product`),
		/// which is evaluated once per group. The contribution of a pair of rows is then the product (PDF row) (coefficients) (FF row)^T.
		/// The observable must be linear in the structure functions, so that it can be computed for each contribution separately.
		/// At NNLO, the points within 1e-8 of the lines xi = 1, zeta = 1, zeta = xi, and zeta = 1 - xi, where the coefficient functions lose precision, contribute nothing.
		/// @tparam pdf_interface_t Type of PDF object. Must satisfy the `is_pdf_interface` concept.
		/// @tparam ff_interface_t Type of FF object. Must satisfy the `is_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @tparam decay_function_t Type of decay function object, used to get the value of the decay function. Must satisfy the `is_decay_function` concept.
		/// @param[in] xi The PDF-side convolution variable
		/// @param[in] zeta The FF-side convolution variable
		/// @param[in] z The FF-side momentum fraction
		/// @param[in] sidis_params Parameters used for the evaluation of the SIDIS structure functions
		/// @param[in] decay_params Parameters used for the evaluation of the decay function
		/// @param[in] F1_qq The quark-to-quark F1 structure function
		/// @param[in] F1_gq The quark-to-gluon F1 structure function
		/// @param[in] F1_qg The gluon-to-quark F1 structure function
		/// @param[in] F1_gg The gluon-to-gluon F1 structure function
		/// @param[in] F2_qq The quark-to-quark F2 structure function
		/// @param[in] F2_gq The quark-to-gluon F2 structure function
		/// @param[in] F2_qg The gluon-to-quark F2 structure function
		/// @param[in] F2_gg The gluon-to-gluon F2 structure function
		/// @param[in] F3_qq The quark-to-quark F3 structure function
		/// @param[in] F3_gq The quark-to-gluon F3 structure function
		/// @param[in] F3_qg The gluon-to-quark F3 structure function
		/// @param[in] F3_gg The gluon-to-gluon F3 structure function
		/// @param[in] observable The observable functor, which gives the observable from the F1, F2, and F3 structure function values
		/// @param[in] fill_component Called as `fill_component(pdf_row, ff_row, physical_channel_index, value)` with each contribution to the observable (see `sidis::component_filler`)
		template <is_pdf_interface pdf_interface_t, is_pdf_interface ff_interface_t, is_coupling coupling_t, is_decay_function decay_function_t>
		void evaluate_observable(
			const double xi, const double zeta, const double z,
			const sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> &sidis_params,
			const decay_parameters<decay_function_t> &decay_params,
			sidis::coefficients::coefficient_function auto &&F1_qq, 
			sidis::coefficients::coefficient_function auto &&F1_gq,
			sidis::coefficients::coefficient_function auto &&F1_qg, 
			sidis::coefficients::coefficient_function auto &&F1_gg,
			sidis::coefficients::coefficient_function auto &&F2_qq, 
			sidis::coefficients::coefficient_function auto &&F2_gq,
			sidis::coefficients::coefficient_function auto &&F2_qg, 
			sidis::coefficients::coefficient_function auto &&F2_gg,
			sidis::coefficients::coefficient_function auto &&F3_qq, 
			sidis::coefficients::coefficient_function auto &&F3_gq,
			sidis::coefficients::coefficient_function auto &&F3_qg, 
			sidis::coefficients::coefficient_function auto &&F3_gg,
			is_observable auto &&observable,
			sidis::component_filler auto &&fill_component
		) const {
			// Checks for thresholds that would lead to NaNs in the coefficients
			if (xi == 1.0 || zeta == 1.0 || xi == zeta) { return; }
			// The NNLO coefficient functions lose precision near the lines xi = 1, zeta = 1, zeta = xi, and zeta = 1 - xi, where their pieces meet or their plus
			// distributions are subtracted: their expressions have poles there that cancel, so that the rounding error of the integrand grows as about eps / d^2
			// with the distance d to a line. At 1e-8 the error is up to about 10% of the integrand, closer to a line the values are meaningless (at 1e-12
			// typically 1e5 times too large), and on it they are not finite. The points within 1e-8 of the lines are left out, and so are those with
			// (1 - xi)(1 - zeta) < 1e-14, since in the corner xi, zeta -> 1 the integrand is not finite below about 3e-16. In the structure functions of the tests,
			// the left-out points make up about 3e-6 of the NNLO value, and up to 4e-5 where it is small from cancellations. At LO and NLO, which don't lose
			// precision, the same points would make up 6e-8 to 5e-7 of the values.
			if (sidis_params.pto == perturbative_order::nnlo) {
				constexpr double line_distance = 1e-8;
				constexpr double corner_product = 1e-14;
				const double one_minus_xi = 1.0 - xi;
				const double one_minus_zeta = 1.0 - zeta;
				if (
					one_minus_xi < line_distance || one_minus_zeta < line_distance || std::abs(zeta - xi) < line_distance
					|| std::abs(one_minus_xi - zeta) < line_distance || one_minus_xi * one_minus_zeta < corner_product
				) { return; }
			}
			// Enforces the cut `zeta < z` for gc coefficients. 
			// For hw coefficients, the check is implemented either explicitly or implicitly directly in the coefficient files.
			if (sidis_params.coefficient_implementation == coefficient_implementation_t::gc && zeta < z) { return; }

			// The PDF and FF values (which depend on the integration point)
			const sidis::parton_values<pdf_interface_t, ff_interface_t, coupling_t> distribution_values{xi, zeta, z, sidis_params, decay_params};
			// A list of values that depend on the integration point but that are shared among the different coefficient functions
			const sidis::evaluation_cache_gc evaluation_cache = make_evaluation_cache(xi, zeta, z, sidis_params);

			// Project out the weight of the F1, F2, and F3 contributions of the observable (which is linear in the Fs)
			const double f1_weight = observable(sidis_params.kinematics, sidis_params.process, 1.0, 0.0, 0.0);
			const double f2_weight = observable(sidis_params.kinematics, sidis_params.process, 0.0, 1.0, 0.0);
			const double f3_weight = observable(sidis_params.kinematics, sidis_params.process, 0.0, 0.0, 1.0);

			for (const auto &[coefficient_channel, group] : coefficient_channels) {
				// The channels of the group have the same coefficient functions, which depend only on their coefficient channel.
				// The coefficients can therefore be evaluated only once per group
				const auto coefficients = [&](auto &&quark_to_quark, auto &&quark_to_gluon, auto &&gluon_to_quark, auto &&gluon_to_gluon) {
					return evaluate_coefficients(
						xi, zeta, z, sidis_params,
						std::forward<decltype(quark_to_quark)>(quark_to_quark),
						std::forward<decltype(quark_to_gluon)>(quark_to_gluon),
						std::forward<decltype(gluon_to_quark)>(gluon_to_quark),
						std::forward<decltype(gluon_to_gluon)>(gluon_to_gluon),
						coefficient_channel, 
						evaluation_cache
					);
				};
				// The coefficients of F1, F2, and F3, evaluated at the first channel with non-zero PDF and FF values. If all of the channels have zero PDF or FF values
				// (e.g. outside the supports of interpolation kernels), the group contributes nothing and the coefficients are never evaluated.
				std::optional<std::array<sidis::vectorization::distribution_product, 3>> structure_function_coefficients;

				for (const partonic_channel_group<sidis::partonic_channel>::indexed_channel &indexed_channel : group) {
					const sidis::partonic_channel &channel = indexed_channel.channel;
					// The PDF and FF values of the channel, which contributes nothing if either of them is zero
					const auto &pdf_rows = distribution_values.pdf_rows(channel, coefficient_channel);
					const auto &ff_rows = distribution_values.ff_rows(channel, coefficient_channel);
					if (std::ranges::empty(pdf_rows) || std::ranges::empty(ff_rows)) { continue; }

					// Actually calculate the coefficients if they haven't been before
					if (!structure_function_coefficients) {
						structure_function_coefficients = std::array{
							coefficients(F1_qq, F1_gq, F1_qg, F1_gg),
							coefficients(F2_qq, F2_gq, F2_qg, F2_gg),
							coefficients(F3_qq, F3_gq, F3_qg, F3_gg)
						};
					}
					const auto &[f1_coefficients, f2_coefficients, f3_coefficients] = *structure_function_coefficients;

					// The coefficients of the observable for this channel, with the couplings and the signs of the antiquark terms in F3
					const double f1_f2_factor = channel_factor(sidis_params, 1, channel);
					const double f3_factor = channel_factor(sidis_params, -1, channel);
					const sidis::vectorization::distribution_product observable_coefficients =
						f1_weight * f1_f2_factor * f1_coefficients + f2_weight * f1_f2_factor * f2_coefficients + f3_weight * f3_factor * f3_coefficients;

					// The contribution of each pair of a row of the PDF values and a row of the FF values: (PDF row) (coefficients) (FF row)^T
					for (const auto &[pdf_row, pdf_values] : pdf_rows) {
						// (PDF row) (coefficients), the coefficient of each FF value, which is the same for every row of the FF values
						const auto ff_coefficients = dm::vectorization::contract(pdf_values, observable_coefficients);
						for (const auto &[ff_row, ff_values] : ff_rows) {
							const double value = dm::vectorization::dot_product(ff_coefficients, ff_values);
							fill_component(pdf_row, ff_row, indexed_channel.physical_channel_index, value);
						}
					}
				}
			}
		}

	private:
		/// @brief Returns the factor that multiplies the value of the coefficient function of \p channel in a structure function:
		/// the coupling, with the sign of the antiquark terms in F3
		/// @param[in] params Parameters used for the evaluation of the SIDIS structure functions
		/// @param[in] sign The "sign" of the structure function, related to the sign of antiquark terms. +1 for F1 and F2, -1 for F3.
		/// @param[in] channel The partonic channel
		static double channel_factor(const auto &params, const int sign, const sidis::partonic_channel &channel) {
			// Compute the value of the coupling, based on the channel
			const double coupling = [&]() {
				const auto &[alpha, beta] = channel.get_coupling_quarks();
				return params.coupling(alpha.flavor, beta.flavor);
			}();

			// Decide the sign, relevant for antiquark channels in F3
			const double coupling_sign = [&]() {
				switch(channel.type) {
					using enum sidis::partonic_channel::channel_type;
					case alphabar_betabar:
					case g_betabar:
					case alphabar_g:
					case alphabarprime_q:
					case betabar_betabar_Fcon:
					case betabar_alpha:
					case beta_betabar_Fcon:
					case qprime_betabar:
					case qbarprime_betabar:
					case alphabarprime_qbar:
						return static_cast<double>(sign);
					default: return 1.0;
				}
			}();

			return coupling_sign * coupling;
		}

		/// Computes the values of the coefficient functions that do not depend on the channel at the given integration point (see `sidis::evaluation_cache_gc`).
		/// Only the gc coefficient functions read them, so for the hw implementation the cache is empty.
		sidis::evaluation_cache_gc make_evaluation_cache(const double xi, const double zeta, const double z, const auto &params) const {
			if (params.coefficient_implementation != coefficient_implementation_t::gc) { return {}; }
			return sidis::evaluation_cache_gc{
				{
					.xi = xi, .zeta = zeta, .z = z,
					.renormalization_scale = params.renormalization_scale,
					.factorization_scale = params.factorization_scale,
					.fragmentation_scale = params.fragmentation_scale,
					.Q2 = params.kinematics.Q2,
					.thresholds = params.thresholds
				}, 
				params.pto
			};
		}

		/// The partonic channels, grouped by their coefficient functions (see `sidis::coefficient_channel`)
		std::flat_map<sidis::coefficient_channel, partonic_channel_group<sidis::partonic_channel>> coefficient_channels;

		/// Returns \p channels grouped by their coefficient functions (see `sidis::coefficient_channel`),
		/// each with the position of its physical channel in \p physical_channels
		static std::flat_map<sidis::coefficient_channel, partonic_channel_group<sidis::partonic_channel>> group_by_coefficients(
			const std::vector<sidis::partonic_channel> &channels,
			const std::flat_set<sidis::physical_channel> &physical_channels
		) {
			std::flat_map<sidis::coefficient_channel, partonic_channel_group<sidis::partonic_channel>> groups;
			for (const sidis::partonic_channel &channel : channels) {
				const auto physical_channel = physical_channels.find(sidis::physical_channel{channel});
				assert(physical_channel != physical_channels.end());
				groups[sidis::coefficient_channel{channel}].push_back({
					channel, 
					static_cast<channel_index>(collections::iterator_to_index(physical_channels, physical_channel))
				});
			}
			return groups;
		}
	};
}

}

#endif