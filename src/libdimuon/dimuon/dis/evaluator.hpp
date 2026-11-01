#ifndef DM_DIS_EVALUATOR_H
#define DM_DIS_EVALUATOR_H

#include <algorithm>
#include <array>
#include <cassert>
#include <flat_map>
#include <flat_set>
#include <optional>
#include <utility>
#include <vector>

#include "dimuon/common/is_observable.hpp"
#include "dimuon/common/partonic_channel_group.hpp"

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/dis/coefficients/coefficient_function.hpp"

#include "dimuon/dis/coefficient_channel.hpp"
#include "dimuon/dis/parameters.hpp"
#include "dimuon/dis/partonic_channel.hpp"
#include "dimuon/dis/parton_values.hpp"
#include "dimuon/dis/physical_channel.hpp"
#include "dimuon/dis/vectorization.hpp"

#include "dimuon/utility/collections.hpp"
#include "dimuon/utility/vectorization.hpp"

namespace dm {

namespace dis {
	/// A component filler is called from `dis::evaluator::evaluate_observable` with
	/// - the physical channel index
	/// - the actual scalar contribution
	/// and is responsible for filling the appropriate integral component with the value
	template <typename function_t>
	concept component_filler = requires(function_t fill_component, const channel_index physical_channel_index, const double value) {
		fill_component(physical_channel_index, value);
	};

	/// Utility class for computing DIS cross section values (see `sidis::evaluator`).
	/// This class is responsible for actually computing integrands.
	class evaluator {
	public:
		/// @brief Constructs the evaluator, groups the partonic channels by their coefficient functions, and indexes each by the position of its physical channel in \p physical_channels
		/// @param[in] channels The partonic channels
		/// @param[in] physical_channels The physical channels (see `dis::physical_channel`), which must include those of \p channels
		explicit evaluator(const std::vector<dis::partonic_channel> &channels, const std::flat_set<dis::physical_channel> &physical_channels) :
			coefficient_channels(group_by_coefficients(channels, physical_channels)) {}

		/// @brief Evaluates the coefficients of the given coefficient function (see `dis::vectorization::distribution_product`)
		/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_scalar_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @param[in] x The Bjorken-x
		/// @param[in] xi The convolution variable
		/// @param[in] params Parameters used for the evaluation of the DIS structure functions
		/// @param[in] integrand The coefficient function of the partonic channels
		/// @param[in] channel The coefficient channel of the partonic channels (see `dis::coefficient_channel`)
		/// @param[in] Q2 The virtuality (GeV^2)
		/// @return The coefficient of each PDF value in the coefficient function (does not include couplings)
		template <is_scalar_pdf_interface pdf_interface_t, is_coupling coupling_t>
		dis::vectorization::distribution_product evaluate_integrand(
			const double x,
			const double xi,
			const dis::parameters<pdf_interface_t, coupling_t> &params,
			dis::coefficients::coefficient_function auto &&integrand,
			const dis::coefficient_channel &channel,
			const double Q2
		) const {
			// Decide the PDF momentum fraction variable based on the mass scheme of the channels
			const double pdf_momentum_fraction = channel.pdf_momentum_fraction(x, Q2);

			// Enforce the lower bound of the xi-integral
			if (xi < pdf_momentum_fraction) { return 0.0; }

			// Compute the coefficients by evaluating the integrand with unit PDF values (see `dis::vectorization::make_masking_vector`)
			using enum dis::vectorization::distribution_value;
			return std::forward<decltype(integrand)>(integrand)({
				.channel = channel,
				.xi = xi,
				.pdf_momentum_fraction = pdf_momentum_fraction,
				.pdf = dis::vectorization::make_masking_vector(pdf),
				.pdf_hat = dis::vectorization::make_masking_vector(pdf_hat),
				.factorization_scale = params.factorization_scale,
				.Q2 = Q2
			});
		}

		/// @brief Evaluates the coefficients of the given structure function for the partonic channels with a given coefficient channel (see `dis::vectorization::distribution_product`)
		/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_scalar_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @param[in] xi The convolution variable
		/// @param[in] params Parameters used for the evaluation of the DIS structure functions
		/// @param[in] quark The quark structure function
		/// @param[in] gluon The gluon structure function
		/// @param[in] channel The coefficient channel of the partonic channels (see `dis::coefficient_channel`)
		/// @return The coefficient of each PDF value (does not include couplings), which is the same for all the channels
		template <is_scalar_pdf_interface pdf_interface_t, is_coupling coupling_t>
		dis::vectorization::distribution_product evaluate_coefficients(
			const double xi,
			const dis::parameters<pdf_interface_t, coupling_t> &params,
			dis::coefficients::coefficient_function auto &&quark,
			dis::coefficients::coefficient_function auto &&gluon,
			const dis::coefficient_channel &channel
		) const {
			// Bound function for evaluating the given integrand
			const auto compute = [&](auto &&integrand) {
				return evaluate_integrand(params.kinematics.x, xi, params, std::forward<decltype(integrand)>(integrand), channel, params.kinematics.Q2);
			};

			// Compute the coefficients, deciding on which provided integrand function to use based on the channel type
			if (channel.is_quark()) {
				return compute(std::forward<decltype(quark)>(quark));
			} else if (channel.is_gluon()) {
				return compute(std::forward<decltype(gluon)>(gluon));
			}
			std::unreachable();
		}

		/// @brief Evaluates the given observable at an integration point, as a sum of contributions: one for each partonic channel
		///
		/// The coefficients of each group of channels are a vector with an element for each PDF value (see `dis::vectorization::distribution_product`),
		/// which is evaluated once per group. The contribution of a channel is then the dot product of its PDF values with the coefficients.
		/// The observable must be linear in the structure functions, so that it can be computed for each contribution separately.
		/// @tparam pdf_interface_t Type of PDF object. Must satisfy the `is_scalar_pdf_interface` concept.
		/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
		/// @param[in] xi The convolution variable
		/// @param[in] dis_params Parameters used for the evaluation of the DIS structure functions
		/// @param[in] F1_q The quark F1 structure function
		/// @param[in] F1_g The gluon F1 structure function
		/// @param[in] F2_q The quark F2 structure function
		/// @param[in] F2_g The gluon F2 structure function
		/// @param[in] F3_q The quark F3 structure function
		/// @param[in] F3_g The gluon F3 structure function
		/// @param[in] observable The observable functor, which gives the observable from the F1, F2, and F3 structure function values
		/// @param[in] fill_component Called as `fill_component(physical_channel_index, value)` with each contribution to the observable (see `dis::component_filler`)
		template <is_scalar_pdf_interface pdf_interface_t, is_coupling coupling_t>
		void evaluate_observable(
			const double xi,
			const dis::parameters<pdf_interface_t, coupling_t> &dis_params,
			dis::coefficients::coefficient_function auto &&F1_q,
			dis::coefficients::coefficient_function auto &&F1_g,
			dis::coefficients::coefficient_function auto &&F2_q,
			dis::coefficients::coefficient_function auto &&F2_g,
			dis::coefficients::coefficient_function auto &&F3_q,
			dis::coefficients::coefficient_function auto &&F3_g,
			is_observable auto &&observable,
			dis::component_filler auto &&fill_component
		) const {
			// Check for the threshold that would lead to NaNs in the coefficients
			if (xi == 1.0) { return; }

			// The PDF values (which depend on the integration point)
			const dis::parton_values<pdf_interface_t, coupling_t> distribution_values{xi, dis_params};

			// Project out the weight of the F1, F2, and F3 contributions of the observable (which is linear in the Fs)
			const double f1_weight = observable(dis_params.kinematics, dis_params.process, 1.0, 0.0, 0.0);
			const double f2_weight = observable(dis_params.kinematics, dis_params.process, 0.0, 1.0, 0.0);
			const double f3_weight = observable(dis_params.kinematics, dis_params.process, 0.0, 0.0, 1.0);

			for (const auto &[coefficient_channel, group] : coefficient_channels) {
				// The channels of the group have the same coefficient functions, which depend only on their coefficient channel.
				// The coefficients can therefore be evaluated only once per group
				const auto coefficients = [&](auto &&quark, auto &&gluon) {
					return evaluate_coefficients(xi, dis_params, std::forward<decltype(quark)>(quark), std::forward<decltype(gluon)>(gluon), coefficient_channel);
				};
				// The coefficients of F1, F2, and F3, evaluated at the first channel with non-zero PDF values. If all of the channels have zero PDF values
				// (e.g. flavors that the PDF doesn't have), the group contributes nothing and the coefficients are never evaluated.
				std::optional<std::array<dis::vectorization::distribution_product, 3>> structure_function_coefficients;

				for (const partonic_channel_group<dis::partonic_channel>::indexed_channel &indexed_channel : group) {
					const dis::partonic_channel &channel = indexed_channel.channel;
					// The PDF values of the channel, which contributes nothing if they are zero
					const std::array<double, dis::vectorization::pdf_value_count> pdf_values = distribution_values.pdf_values(channel, coefficient_channel);
					if (std::ranges::all_of(pdf_values, [](const double value) { return value == 0.0; })) { continue; }

					// Actually calculate the coefficients if they haven't been before
					if (!structure_function_coefficients) {
						structure_function_coefficients = std::array{
							coefficients(F1_q, F1_g),
							coefficients(F2_q, F2_g),
							coefficients(F3_q, F3_g)
						};
					}
					const auto &[f1_coefficients, f2_coefficients, f3_coefficients] = *structure_function_coefficients;

					// The coefficients of the observable for this channel, with the couplings and the signs of the antiquark terms in F3
					const double f1_f2_factor = channel_factor(dis_params, 1, channel);
					const double f3_factor = channel_factor(dis_params, -1, channel);
					const dis::vectorization::distribution_product observable_coefficients =
						f1_weight * f1_f2_factor * f1_coefficients + f2_weight * f1_f2_factor * f2_coefficients + f3_weight * f3_factor * f3_coefficients;

					fill_component(indexed_channel.physical_channel_index, dm::vectorization::dot_product(pdf_values, observable_coefficients));
				}
			}
		}

	private:
		/// @brief Returns the factor that multiplies the value of the coefficient function of \p channel in a structure function:
		/// the coupling, with the sign of the antiquark terms in F3
		/// @param[in] params Parameters used for the evaluation of the DIS structure functions
		/// @param[in] sign The "sign" of the structure function, related to the sign of antiquark terms. +1 for F1 and F2, -1 for F3.
		/// @param[in] channel The partonic channel
		static double channel_factor(const auto &params, const int sign, const dis::partonic_channel &channel) {
			// Compute the value of the coupling, based on the channel
			const double coupling = [&]() {
				const auto &[alpha, beta] = channel.get_coupling_quarks();
				return params.coupling(alpha.flavor, beta.flavor);
			}();

			// Decide the sign, relevant for the charge-conjugate channels in F3
			const double coupling_sign = [&]() {
				switch (channel.type) {
					using enum dis::partonic_channel::channel_type;
					case alphabar_betabar:
					case alphabar_g:
						return static_cast<double>(sign);
					default: return 1.0;
				}
			}();

			return coupling_sign * coupling;
		}

		/// The partonic channels, grouped by their coefficient functions (see `dis::coefficient_channel`)
		std::flat_map<dis::coefficient_channel, partonic_channel_group<dis::partonic_channel>> coefficient_channels;

		/// Returns \p channels grouped by their coefficient functions (see `dis::coefficient_channel`),
		/// each with the position of its physical channel in \p physical_channels
		static std::flat_map<dis::coefficient_channel, partonic_channel_group<dis::partonic_channel>> group_by_coefficients(
			const std::vector<dis::partonic_channel> &channels,
			const std::flat_set<dis::physical_channel> &physical_channels
		) {
			std::flat_map<dis::coefficient_channel, partonic_channel_group<dis::partonic_channel>> groups;
			for (const dis::partonic_channel &channel : channels) {
				const auto physical_channel = physical_channels.find(dis::physical_channel{channel});
				assert(physical_channel != physical_channels.end());
				groups[dis::coefficient_channel{channel}].push_back({
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
