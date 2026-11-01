#ifndef DM_SIDIS_INSTANCE_H
#define DM_SIDIS_INSTANCE_H

#include <cstddef>
#include <flat_set>
#include <format>
#include <ostream>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "dimuon/common/coefficient_implementation.hpp"
#include "dimuon/common/heavy_quark_thresholds.hpp"
#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/common/perturbative_quantity.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/common/is_observable.hpp"
#include "dimuon/common/cross_section_functions.hpp"

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/decay/decay_instance.hpp"

#include "dimuon/sidis/decay_parameters.hpp"
#include "dimuon/sidis/evaluator.hpp"
#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/sidis/physical_channel.hpp"
#include "dimuon/sidis/parameters.hpp"
#include "dimuon/sidis/result.hpp"

#include "dimuon/integration/integrator_mc.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/utility/shims.hpp"
#include "dimuon/utility/tensor.hpp"

#include "dimuon/common/scale_dependence.hpp"

#include "dimuon/sidis/coefficients/f1_hw.hpp"
#include "dimuon/sidis/coefficients/f2_hw.hpp"
#include "dimuon/sidis/coefficients/f3_hw.hpp"
#include "dimuon/sidis/coefficients/f1_gc.hpp"
#include "dimuon/sidis/coefficients/f2_gc.hpp"
#include "dimuon/sidis/coefficients/f3_gc.hpp"

namespace dm {

namespace sidis {
	/// @brief Entry point for computing SIDIS cross section values
	/// @details The configuration is checked once, at construction, which throws `std::invalid_argument` if it is invalid (see `check_configuration`).
	/// @tparam pdf_interface_t Type of PDF object, used for obtaining PDF values. Must satisfy the `is_pdf_interface` concept.
	/// @tparam ff_interface_t Type of FF object, used for obtaining FF values. Must satisfy the `is_pdf_interface` concept.
	/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`), used to get the coupling between partons. Must satisfy the `is_coupling` concept.
	/// @tparam decay_function_t Type of decay function object, used to get the value of the decay function. Must satisfy the `is_decay_function` concept.
	/// @tparam renormalization_scale_t Type of renormalization scale function object. Must satisfy the `is_scale_dependence` concept.
	/// @tparam factorization_scale_t Type of factorization scale function object. Must satisfy the `is_scale_dependence` concept.
	/// @tparam fragmentation_scale_t Type of fragmentation scale function object. Must satisfy the `is_scale_dependence` concept.
	/// @tparam observable_t Type of observable functor. Must satisfy the `is_observable` concept.
	template <
		is_pdf_interface pdf_interface_t, 
		is_pdf_interface ff_interface_t, 
		is_coupling coupling_t,
		is_decay_function decay_function_t,
		is_scale_dependence renormalization_scale_t,
		is_scale_dependence factorization_scale_t, 
		is_scale_dependence fragmentation_scale_t,
		is_observable observable_t
	>
	class instance {
	public:
		/// Coupling object, used to get the coupling between partons 
		coupling_t coupling;
		/// PDF object, used for obtaining PDF values
		pdf_interface_t pdf;
		/// FF objects, used for obtaining FF values: the FF of each decay channel, in the same order as `decays`
		std::vector<ff_interface_t> ffs;
		/// The decay (resonance and decay function) of each decay channel, in the same order as `ffs`
		std::vector<decay_instance<decay_function_t>> decays;
		/// Integration parameters used for all integrals
		std::unordered_map<perturbative_order, dm::integration_parameters> integration_parameters;
		/// The scattering process
		scattering_process process;
		/// Function object for the renormalization scale
		scale_dependence::function<renormalization_scale_t> renormalization_scale_function;
		/// Function object for the factorization scale
		scale_dependence::function<factorization_scale_t> factorization_scale_function;
		/// Function object for the fragmentation scale
		scale_dependence::function<fragmentation_scale_t> fragmentation_scale_function;
		/// Flag for the cross section prefactor - see `cross_section_functions::prefactor`
		bool use_modified_cross_section_prefactor;
		/// Function object that defines the actual observable being computed from structure functions
		observable_t observable;
		/// Highest perturbative order to be calculated (i.e. all orders up to and including this are calculated)
		perturbative_order order;
		/// Minimum energy of the primary muon (GeV). Determines the maximum allowed value of y.
		double primary_muon_min_energy;
		/// List of partonic channels per perturbative order
		const std::unordered_map<perturbative_order, std::vector<sidis::partonic_channel>> channels;
		/// Coefficient implementation per perturbative order, see `coefficient_implementation_t`.
		/// The LO and NLO coefficient functions have both implementations, and the NNLO ones only the `gc` implementation.
		std::unordered_map<perturbative_order, coefficient_implementation_t> coefficient_implementation;
		/// Thresholds of the heavy quarks, which determine the number of active flavors
		/// in the NNLO coefficient functions. Usually the quark masses of the PDF set in use.
		heavy_quark_thresholds thresholds;
		/// Output stream where progress reporting is printed. Defaults to `std::cout`. 
		/// Progress reporting can be disabled by setting this to `io::null_stream`.
		std::ostream &progress_stream = std::cout;

		/// Set of physical channels; see the section "Channels" in `numerical_evaluation.md` for an explanation.
		/// Its initializer first checks the configuration, once, at construction (see `check_configuration`).
		/// The check reads the other members, so this must stay the last member.
		const std::flat_set<sidis::physical_channel> physical_channels = [this] {
			check_configuration();
			return channels
				| std::views::filter([this](const auto &order_channels) { return order_channels.first <= order; })
				| std::views::values
				| std::views::join
				| std::views::transform([](const sidis::partonic_channel &channel) { return sidis::physical_channel{channel}; })
				| std::ranges::to<std::flat_set>();
		}();

	private:
		/// Checks that the configuration, with which this instance was constructed from, is valid
		void check_configuration() const {
			for (const perturbative_order calculation_order : {perturbative_order::lo, perturbative_order::nlo, perturbative_order::nnlo}) {
				if (calculation_order > order) { break; }
				if (!channels.contains(calculation_order)) {
					throw std::invalid_argument(std::format("The channels of {} are missing: each order up to {} needs a list of channels, which may be empty", calculation_order, order));
				}
				if (!coefficient_implementation.contains(calculation_order)) {
					throw std::invalid_argument(std::format("The coefficient implementation of {} is missing: each order up to {} needs one", calculation_order, order));
				}
				if (!integration_parameters.contains(calculation_order)) {
					throw std::invalid_argument(std::format("The integration parameters of {} are missing: each order up to {} needs them", calculation_order, order));
				}
			}

			const auto nnlo_implementation = coefficient_implementation.find(perturbative_order::nnlo);
			if (nnlo_implementation != coefficient_implementation.end() && nnlo_implementation->second == coefficient_implementation_t::hw) {
				throw std::invalid_argument("The hw coefficients are not available at NNLO");
			}

			for (const auto &[channels_order, order_channels] : channels) {
				// The mass terms of the SACOT schemes are only in the hw coefficient functions of NLO (at LO, the SACOT schemes are those of slow rescaling)
				const auto implementation = coefficient_implementation.find(channels_order);
				const bool gc_at_nlo = channels_order == perturbative_order::nlo
					&& implementation != coefficient_implementation.end() && implementation->second == coefficient_implementation_t::gc;

				for (const sidis::partonic_channel &channel : order_channels) {
					if (channel.pto != channels_order) {
						throw std::invalid_argument(std::format("A channel of {} is listed with the channels of {}: {}", channel.pto, channels_order, channel));
					}
					if (!channel.check_charge_conservation(process)) {
						throw std::invalid_argument(std::format("A channel doesn't conserve the electric charge in the process: {}", channel));
					}
					if (is_sacot_scheme(channel.mass_scheme) && channels_order == perturbative_order::nnlo) {
						throw std::invalid_argument(std::format("The SACOT schemes are not available at NNLO: {}", channel));
					}
					if (is_sacot_scheme(channel.mass_scheme) && gc_at_nlo) {
						throw std::invalid_argument(std::format("The SACOT schemes are not available with the gc coefficients at NLO: {}", channel));
					}
				}
			}

			if (decays.empty() || ffs.size() != decays.size()) {
				throw std::invalid_argument(std::format(
					"Each decay channel needs an FF, and there must be at least one decay channel (decay channels: {}, FFs: {})", decays.size(), ffs.size()
				));
			}
		}

		/// Calls \p evaluate with the appropriate coefficient functions,
		/// which are selected by the perturbative order \p calculation_order and the associated coefficient implementation (see `coefficient_implementation`)
		decltype(auto) with_coefficient_functions(const perturbative_order calculation_order, auto &&evaluate) const {
			using namespace sidis::coefficients;

			const coefficient_implementation_t impl = coefficient_implementation.at(calculation_order);
			switch(calculation_order) {
				case perturbative_order::lo:
				{
					switch(impl) {
						case coefficient_implementation_t::hw: 
							return evaluate(
								f1_hw::lo::quark_to_quark, f1_hw::lo::quark_to_gluon, f1_hw::lo::gluon_to_quark, f1_hw::lo::gluon_to_gluon,
								f2_hw::lo::quark_to_quark, f2_hw::lo::quark_to_gluon, f2_hw::lo::gluon_to_quark, f2_hw::lo::gluon_to_gluon,
								f3_hw::lo::quark_to_quark, f3_hw::lo::quark_to_gluon, f3_hw::lo::gluon_to_quark, f3_hw::lo::gluon_to_gluon
							);
						case coefficient_implementation_t::gc: 
							return evaluate(
								f1_gc::lo::quark_to_quark, f1_gc::lo::quark_to_gluon, f1_gc::lo::gluon_to_quark, f1_gc::lo::gluon_to_gluon,
								f2_gc::lo::quark_to_quark, f2_gc::lo::quark_to_gluon, f2_gc::lo::gluon_to_quark, f2_gc::lo::gluon_to_gluon,
								f3_gc::lo::quark_to_quark, f3_gc::lo::quark_to_gluon, f3_gc::lo::gluon_to_quark, f3_gc::lo::gluon_to_gluon
							);
					}
					std::unreachable();
				}
				case perturbative_order::nlo:
				{
					switch(impl) {
						case coefficient_implementation_t::hw: 
							return evaluate(
								f1_hw::nlo::quark_to_quark, f1_hw::nlo::quark_to_gluon, f1_hw::nlo::gluon_to_quark, f1_hw::nlo::gluon_to_gluon, 
								f2_hw::nlo::quark_to_quark, f2_hw::nlo::quark_to_gluon, f2_hw::nlo::gluon_to_quark, f2_hw::nlo::gluon_to_gluon, 
								f3_hw::nlo::quark_to_quark, f3_hw::nlo::quark_to_gluon, f3_hw::nlo::gluon_to_quark, f3_hw::nlo::gluon_to_gluon
							);
						case coefficient_implementation_t::gc: 
							return evaluate(
								f1_gc::nlo::quark_to_quark, f1_gc::nlo::quark_to_gluon, f1_gc::nlo::gluon_to_quark, f1_gc::nlo::gluon_to_gluon,
								f2_gc::nlo::quark_to_quark, f2_gc::nlo::quark_to_gluon, f2_gc::nlo::gluon_to_quark, f2_gc::nlo::gluon_to_gluon,
								f3_gc::nlo::quark_to_quark, f3_gc::nlo::quark_to_gluon, f3_gc::nlo::gluon_to_quark, f3_gc::nlo::gluon_to_gluon
							);
					}
					std::unreachable();
				}
				case perturbative_order::nnlo:
				{
					switch(impl) {
						case coefficient_implementation_t::hw:
							// There are no hw NNLO coefficient functions, which `check_configuration` rules out
							std::unreachable();
						case coefficient_implementation_t::gc:
							return evaluate(
								f1_gc::nnlo::quark_to_quark, f1_gc::nnlo::quark_to_gluon, f1_gc::nnlo::gluon_to_quark, f1_gc::nnlo::gluon_to_gluon,
								f2_gc::nnlo::quark_to_quark, f2_gc::nnlo::quark_to_gluon, f2_gc::nnlo::gluon_to_quark, f2_gc::nnlo::gluon_to_gluon,
								f3_gc::nnlo::quark_to_quark, f3_gc::nnlo::quark_to_gluon, f3_gc::nnlo::gluon_to_quark, f3_gc::nnlo::gluon_to_gluon
							);
					}
					std::unreachable();
				}
			}
			std::unreachable();
		}

		/// Returns the parameters for evaluating the SIDIS structure functions of \p calculation_order at the kinematical point \p kinematics
		sidis::parameters<pdf_interface_t, ff_interface_t, coupling_t> make_sidis_parameters(const trf_kinematics &kinematics, const perturbative_order calculation_order) const {
			return sidis::parameters{
				pdf,
				ffs,
				coupling, process, kinematics,
				calculation_order,
				channels.at(calculation_order),
				coefficient_implementation.at(calculation_order),
				thresholds,
				renormalization_scale_function(kinematics), 
				factorization_scale_function(kinematics), 
				fragmentation_scale_function(kinematics)
			};
		}

	public:
		/// Computes the minimum `z_min` across all decay channels at the given kinematical point
		double compute_global_z_min(const trf_kinematics &kinematics) const noexcept {
			return std::ranges::min(decays | std::views::transform([&kinematics](const auto &decay) { return decay.z_min(kinematics); }));
		}

		/// Computes the value of `alpha_s` at the given kinematical point from the PDF (zero for the interpolation kernels, see `vector_kernel_interface::alpha_s`)
		double compute_alpha_s(const trf_kinematics &kinematics) const {
			const double renormalization_scale = renormalization_scale_function(kinematics);
			return pdf.alpha_s(renormalization_scale);
		}

		/// Computes the dimuon cross section `d^2 sigma/dx dQ^2` at the given kinematical point
		sidis::result lepton_pair_cross_section_xQ2(const trf_kinematics &kinematics) const {
			const auto integrate_order = [&](
				const perturbative_order calculation_order, const auto &integrand, const std::vector<double> &lower, const std::vector<double> &upper, const std::size_t n_components
			) {
				const integrator_mc integrator(integration_parameters.at(calculation_order));
				return integrator.integrate(integrand, lower, upper, n_components);
			};

			return integrate_common(kinematics, integrate_order);
		}

		/// Computes the dimuon cross section `d^2 sigma/dx dy` at the given kinematical point
		sidis::result lepton_pair_cross_section_xy(const trf_kinematics &kinematics) const {
			const double jacobian = cross_section_functions::xy_jacobian(kinematics, process);
			return lepton_pair_cross_section_xQ2(kinematics) * jacobian;
		}

	private:
		/// @brief Integrates the cross section `d^2 sigma/dx dQ^2` at each perturbative order up to `order`, and collects the results.
		/// The components of the integral of each order are the elements of the channel decomposition: a tensor over the rows of the PDF values,
		/// the rows of the FF values, and the physical channels (see `sidis::evaluator::evaluate_observable` and `physical_channels`).
		/// @param[in] kinematics The kinematical point
		/// @param[in] integrate_order Integrates one perturbative order, called as `integrate_order(calculation_order, integrand, lower, upper, n_components)`:
		/// `integrand` is the integrand of the cross section, called as `integrand(input, components)` at the integration point `input`,
		/// which adds each contribution to its component, `components[component] += value` (see `is_mc_integrand_function`),
		/// `lower` and `upper` are the bounds of the integral, and `n_components` is the number of components. Returns an `integrator_mc::result`.
		/// @return The totals and the channel decomposition, multiplied by the cross section prefactor
		sidis::result integrate_common(const trf_kinematics &kinematics, auto &&integrate_order) const {
			const double z_min = compute_global_z_min(kinematics);
			// Sanity check z_min
			if (z_min < 0.0 || z_min > 1.0) { return {}; }

			const std::size_t pdf_row_count = pdf.dimension();
			const std::size_t ff_row_count = std::ranges::max(ffs | std::views::transform([](const ff_interface_t &ff) { return ff.dimension(); }));

			// Layout of the channel decomposition, which is also the layout of the components of the integral of each perturbative order
			const shims::layout_right::mapping layout{shims::dextents<std::size_t, 3>{pdf_row_count, ff_row_count, physical_channels.size()}};

			// The integral of each perturbative order sets the values of that order
			const double alpha_s = compute_alpha_s(kinematics);
			sidis::result result{perturbative_quantity{alpha_s}, physical_channels, tensor<perturbative_quantity, 3>{layout.extents(), perturbative_quantity{alpha_s}}};

			// Integrate the cross section at the requested perturbative order
			const auto perform_integral = [&](const perturbative_order calculation_order) {
				// If there are no channels, the values stay zero
				const std::vector<sidis::partonic_channel> &order_channels = channels.at(calculation_order);
				if (order_channels.empty()) { return; }

				// Enforce phase space cuts for y (see `primary_muon_min_energy`): above the cut, the values stay zero
				const double y_max = 1.0 - primary_muon_min_energy / kinematics.E_beam;
				if (kinematics.y > y_max) { return; }

				// Construct necessary parameter objects
				const sidis::parameters sidis_params = make_sidis_parameters(kinematics, calculation_order);
				const decay_parameters decay_params{ decays };
				const sidis::evaluator evaluator{order_channels, physical_channels};

				// The integrand as a function of the integration point, which is (xi, zeta, z): the contributions to the cross section
				// (see `sidis::evaluator::evaluate_observable`), each added to its component, in the layout of the channel decomposition
				const auto integrand = [&](const std::span<const double> input, const std::span<double> components) {
					// Evaluate with the coefficient functions of the perturbative order and coefficient implementation
					with_coefficient_functions(calculation_order, [&](auto &&...coefficient_functions) {
						evaluator.evaluate_observable(
							input[0], input[1], input[2], sidis_params, decay_params,
							std::forward<decltype(coefficient_functions)>(coefficient_functions)...,
							observable,
							[&](
								const distribution_index pdf_row,
								const distribution_index ff_row,
								const channel_index physical_channel_index,
								const double value
							) {
								components[layout(pdf_row, ff_row, physical_channel_index)] += value;
							}
						);
					});
				};
				const integrator_mc::result integral = integrate_order(
					calculation_order, integrand, {kinematics.x, 0.0, z_min}, {1.0, 1.0, 1.0}, layout.required_span_size()
				);

				result.value.raw_correction(calculation_order) = integral.value();
				// The components are the elements of the channel decomposition, in the same (row-major) order
				for (auto &&[element, component] : std::views::zip(result.channel_decomposition, integral.component_values())) {
					element.raw_correction(calculation_order) = component;
				}
			};

			// Integrates the order `calculation_order` and reports the progress: `message` before the integral and "done" after it.
			// The progress stream is flushed after each, so that the progress shows while the integral is computed.
			const auto perform_integral_with_progress = [&](const perturbative_order calculation_order, const std::string_view message) {
				std::print(progress_stream, "{}", message);
				progress_stream.flush();

				perform_integral(calculation_order);

				std::println(progress_stream, "done");
				progress_stream.flush();
			};

			perform_integral_with_progress(perturbative_order::lo, "Calculating LO...   ");
			if (order >= perturbative_order::nlo) {
				perform_integral_with_progress(perturbative_order::nlo, "Calculating NLO...  ");
			}
			if (order >= perturbative_order::nnlo) {
				perform_integral_with_progress(perturbative_order::nnlo, "Calculating NNLO... ");
			}

			// Multiply result (including channel decompositions) by the cross section prefactor
			const double prefactor = cross_section_functions::prefactor(kinematics, use_modified_cross_section_prefactor);

			return prefactor * std::move(result);
		}
	};
}

}

#endif