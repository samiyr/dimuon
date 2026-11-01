#ifndef DM_DIS_INSTANCE_H
#define DM_DIS_INSTANCE_H

#include <flat_set>
#include <format>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "dimuon/common/cross_section_functions.hpp"
#include "dimuon/common/is_observable.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/perturbative_quantity.hpp"
#include "dimuon/common/scale_dependence.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/structure_function.hpp"
#include "dimuon/common/trf_kinematics.hpp"

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/integration/integration_parameters.hpp"
#include "dimuon/integration/integrator_cquad.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/dis/evaluator.hpp"
#include "dimuon/dis/parameters.hpp"
#include "dimuon/dis/partonic_channel.hpp"
#include "dimuon/dis/physical_channel.hpp"
#include "dimuon/dis/result.hpp"

#include "dimuon/dis/coefficients/f1.hpp"
#include "dimuon/dis/coefficients/f2.hpp"
#include "dimuon/dis/coefficients/f3.hpp"

namespace dm {

namespace dis {
	/// @brief Entry point for computing DIS cross section and structure function values (see `sidis::instance`)
	/// @details The configuration is checked once, at construction, which throws `std::invalid_argument` if it is invalid (see `check_configuration`).
	/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_scalar_pdf_interface` concept.
	/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
	/// @tparam renormalization_scale_t Type of renormalization scale function object. Must satisfy the `is_scale_dependence` concept.
	/// @tparam factorization_scale_t Type of factorization scale function object. Must satisfy the `is_scale_dependence` concept.
	/// @tparam observable_t Type of observable functor. Must satisfy the `is_observable` concept.
	template <
		is_scalar_pdf_interface pdf_interface_t,
		is_coupling coupling_t,
		is_scale_dependence renormalization_scale_t,
		is_scale_dependence factorization_scale_t,
		is_observable observable_t
	>
	class instance {
	public:
		/// Coupling object, used to get the coupling between partons
		coupling_t coupling;
		/// PDF object, used for obtaining PDF values
		pdf_interface_t pdf;
		/// Integration parameters per perturbative order, of which the CQUAD parameters are used
		std::unordered_map<perturbative_order, dm::integration_parameters> integration_parameters;
		/// The scattering process
		scattering_process process;
		/// Function object for the renormalization scale
		scale_dependence::function<renormalization_scale_t> renormalization_scale_function;
		/// Function object for the factorization scale
		scale_dependence::function<factorization_scale_t> factorization_scale_function;
		/// Flag for the cross section prefactor - see `cross_section_functions::prefactor`
		bool use_modified_cross_section_prefactor;
		/// Function object that defines the observable of the cross section, computed from structure functions
		observable_t observable;
		/// Highest perturbative order to be calculated (i.e. all orders up to and including this are calculated). There are DIS coefficient functions up to NLO.
		perturbative_order order;
		/// Minimum energy of the primary muon (GeV). Determines the maximum allowed value of y.
		double primary_muon_min_energy;
		/// List of partonic channels per perturbative order
		const std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> channels;

		/// Set of physical channels, the initial-state flavors of the channels (see `dis::physical_channel`).
		/// Its initializer first checks the configuration, once, at construction (see `check_configuration`).
		/// The check reads the other members, so this must stay the last member.
		const std::flat_set<dis::physical_channel> physical_channels = [this] {
			check_configuration();
			return channels
				| std::views::filter([this](const auto &order_channels) { return order_channels.first <= order; })
				| std::views::values
				| std::views::join
				| std::views::transform([](const dis::partonic_channel &channel) { return dis::physical_channel{channel}; })
				| std::ranges::to<std::flat_set>();
		}();

	private:
		/// Checks that the configuration, with which this instance was constructed from, is valid
		void check_configuration() const {
			if (order > perturbative_order::nlo) {
				throw std::invalid_argument("There are DIS coefficient functions only up to NLO");
			}

			for (const perturbative_order calculation_order : {perturbative_order::lo, perturbative_order::nlo}) {
				if (calculation_order > order) { break; }
				if (!channels.contains(calculation_order)) {
					throw std::invalid_argument(std::format("The channels of {} are missing: each order up to {} needs a list of channels, which may be empty", calculation_order, order));
				}
				if (!integration_parameters.contains(calculation_order)) {
					throw std::invalid_argument(std::format("The integration parameters of {} are missing: each order up to {} needs them", calculation_order, order));
				}
			}

			for (const auto &[channels_order, order_channels] : channels) {
				for (const dis::partonic_channel &channel : order_channels) {
					if (channel.pto != channels_order) {
						throw std::invalid_argument(std::format("A channel of {} is listed with the channels of {}: {}", channel.pto, channels_order, channel));
					}
					if (!channel.check_charge_conservation(process)) {
						throw std::invalid_argument(std::format("A channel doesn't conserve the electric charge in the process: {}", channel));
					}
				}
			}
		}

		/// Calls \p evaluate with the coefficient functions of the perturbative order \p calculation_order (see `dis::evaluator::evaluate_observable`)
		static decltype(auto) with_coefficient_functions(const perturbative_order calculation_order, auto &&evaluate) {
			using namespace dis::coefficients;

			switch (calculation_order) {
				case perturbative_order::lo:
					return evaluate(f1::lo::quark, f1::lo::gluon, f2::lo::quark, f2::lo::gluon, f3::lo::quark, f3::lo::gluon);
				case perturbative_order::nlo:
					return evaluate(f1::nlo::quark, f1::nlo::gluon, f2::nlo::quark, f2::nlo::gluon, f3::nlo::quark, f3::nlo::gluon);
				case perturbative_order::nnlo:
					// There are no NNLO DIS coefficient functions, which `check_configuration` rules out
					std::unreachable();
			}
			std::unreachable();
		}

		/// Returns the parameters for evaluating the DIS structure functions of the partonic channels \p order_channels at the kinematical point \p kinematics
		dis::parameters<pdf_interface_t, coupling_t> make_dis_parameters(const trf_kinematics &kinematics, const std::vector<dis::partonic_channel> &order_channels) const {
			return dis::parameters{pdf, coupling, process, kinematics, order_channels, factorization_scale_function(kinematics)};
		}

	public:
		/// Computes the value of `alpha_s` at the given kinematical point from the PDF
		double compute_alpha_s(const trf_kinematics &kinematics) const {
			const double renormalization_scale = renormalization_scale_function(kinematics);
			return pdf.alpha_s(renormalization_scale);
		}

		/// Computes the cross section `d^2 sigma/dx dQ^2` of the observable `observable` at the given kinematical point
		dis::result differential_cross_section_xQ2(const trf_kinematics &kinematics) const {
			// Enforce phase space cuts for y (see `primary_muon_min_energy`): above the cut, the values stay zero
			const double y_max = 1.0 - primary_muon_min_energy / kinematics.E_beam;
			if (kinematics.y > y_max) { return make_result(kinematics); }

			const double prefactor = cross_section_functions::prefactor(kinematics, use_modified_cross_section_prefactor);
			return prefactor * integrate_common(kinematics, observable);
		}

		/// Computes the cross section `d^2 sigma/dx dy` of the observable `observable` at the given kinematical point
		dis::result differential_cross_section_xy(const trf_kinematics &kinematics) const {
			const double jacobian = cross_section_functions::xy_jacobian(kinematics, process);
			return differential_cross_section_xQ2(kinematics) * jacobian;
		}

		/// Computes the structure function \p sf at the given kinematical point
		template <structure_function_t sf>
		dis::result structure_function(const trf_kinematics &kinematics) const {
			return integrate_common(kinematics, structure_function::observable<sf>{});
		}

	private:
		/// Returns a result with zero values at the kinematical point \p kinematics, with an element of the channel decomposition for each physical channel
		dis::result make_result(const trf_kinematics &kinematics) const {
			const double alpha_s = compute_alpha_s(kinematics);
			return dis::result{perturbative_quantity{alpha_s}, physical_channels, std::vector(physical_channels.size(), perturbative_quantity{alpha_s})};
		}

		/// @brief Integrates \p observable_ (a cross section or a structure function) at each perturbative order up to `order`, and collects the results.
		/// Each physical channel is integrated on its own with CQUAD, which integrates one function at a time: the integral of the observable over the convolution variable xi
		/// of the channels of the physical channel (see `dis::evaluator::evaluate_observable`). The integrand of a channel starts at its momentum fraction, so the range of xi is split
		/// at the momentum fractions. The total is the sum of the channel decomposition.
		/// @param[in] kinematics The kinematical point
		/// @param[in] observable_ The observable, which gives the integrand from the structure functions
		/// @return The totals and the channel decomposition
		dis::result integrate_common(const trf_kinematics &kinematics, const is_observable auto &observable_) const {
			dis::result result = make_result(kinematics);

			// Integrate the observable at the requested perturbative order
			const auto perform_integral = [&](const perturbative_order calculation_order) {
				const std::vector<dis::partonic_channel> &order_channels = channels.at(calculation_order);
				const integrator_cquad integrator(integration_parameters.at(calculation_order));

				double total = 0.0;
				for (auto &&[component, element] : std::views::zip(physical_channels, result.channel_decomposition)) {
					// The partonic channels of the physical channel `component`. If there are none, the value stays zero
					const std::vector<dis::partonic_channel> component_channels = order_channels
						| std::views::filter([&](const dis::partonic_channel &channel) { return dis::physical_channel{channel} == component; })
						| std::ranges::to<std::vector>();
					if (component_channels.empty()) { continue; }

					// Construct necessary parameter objects
					const dis::parameters dis_params = make_dis_parameters(kinematics, component_channels);
					const dis::evaluator evaluator{component_channels, physical_channels};

					// The integrand as a function of the convolution variable xi: the sum of the contributions of the channels (see `dis::evaluator::evaluate_observable`)
					const auto integrand = [&](const double xi) {
						double value = 0.0;
						with_coefficient_functions(calculation_order, [&](auto &&...coefficient_functions) {
							evaluator.evaluate_observable(
								xi, dis_params,
								std::forward<decltype(coefficient_functions)>(coefficient_functions)...,
								observable_,
								[&](const channel_index, const double contribution) { value += contribution; }
							);
						});
						return value;
					};

					// The integral from the smallest momentum fraction to one, split at the momentum fractions (those at or above one are outside the range)
					std::vector<double> edges = dis_params.momentum_fractions
						| std::views::filter([](const double momentum_fraction) { return momentum_fraction < 1.0; })
						| std::ranges::to<std::vector>();
					edges.push_back(1.0);

					double integral = 0.0;
					for (const auto &[lower, upper] : edges | std::views::pairwise) {
						integral += integrator.integrate(integrand, lower, upper).value;
					}

					element.raw_correction(calculation_order) = integral;
					total += integral;
				}

				result.value.raw_correction(calculation_order) = total;
			};

			perform_integral(perturbative_order::lo);
			if (order >= perturbative_order::nlo) {
				perform_integral(perturbative_order::nlo);
			}

			return result;
		}
	};
}

}

#endif
