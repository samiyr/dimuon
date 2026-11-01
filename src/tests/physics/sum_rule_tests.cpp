#include <array>
#include <cmath>
#include <flat_set>
#include <format>
#include <ranges>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/structure_function.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/dis/channel_factory.hpp>
#include <dimuon/dis/coefficients/f1.hpp>
#include <dimuon/dis/coefficients/f2.hpp>
#include <dimuon/dis/coefficients/f3.hpp>
#include <dimuon/dis/evaluator.hpp>
#include <dimuon/dis/parameters.hpp>
#include <dimuon/dis/partonic_channel.hpp>
#include <dimuon/dis/physical_channel.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/integration/integrator_mc.hpp>
#include <dimuon/pdf/interface/functional_form_interface.hpp>
#include <dimuon/sidis/channel_factory.hpp>
#include <dimuon/sidis/evaluator.hpp>

#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum structure_function_t;

	const double x = 0.1;
	const double Q2 = 10.0;
	const scattering_process process{scattering_process::process_type::neutrino, constants::particles::proton, constants::particles::neutrino};
	const parton_configuration configuration{std::vector{flavor::u, flavor::d, flavor::s, flavor::c}, std::array<double, constants::total_flavors>{}};
	const channel_configuration zero_mass{mass_scheme_t::zero_mass};

	template <structure_function_t sf>
	double dis_integrand(const double xi, const perturbative_order pto, const dis::evaluator &evaluator, const auto &params) {
		using namespace dis::coefficients;
		double total = 0.0;
		const auto add = [&](const channel_index, const double value) { total += value; };
		if (pto == lo) {
			evaluator.evaluate_observable(xi, params, f1::lo::quark, f1::lo::gluon, f2::lo::quark, f2::lo::gluon, f3::lo::quark, f3::lo::gluon, structure_function::observable<sf>{}, add);
		} else {
			evaluator.evaluate_observable(xi, params, f1::nlo::quark, f1::nlo::gluon, f2::nlo::quark, f2::nlo::gluon, f3::nlo::quark, f3::nlo::gluon, structure_function::observable<sf>{}, add);
		}
		return total;
	}

	/// SIDIS integrates to DIS
	void expect_momentum_sum_rule(const perturbative_order pto) {
		const trf_kinematics kinematics = trf_kinematics::Q2_s(x, Q2, 1000.0, process.target.mass, process.projectile.mass);
		const ckm_matrix coupling = ckm_matrix::global_fit();
		const auto pdf = toy_pdf_interface();
		// A single "hadron" with all the momentum of every parton: z D(z) = 6 z (1 - z) integrates to one
		const functional_form_interface ff{[](const flavor::type, const double z, const double) { return 6.0 * z * (1.0 - z); }};
		const auto weight_z = [](const double, const double z, const double, const particle &) { return z; };
		const decay_parameters decays{std::vector{decay_instance{particle{0.0, 1.0}, constants::particles::proton, weight_z, 0.0}}};

		const std::vector<dis::partonic_channel> dis_channels = dis::channel_factory::create_all_channels(pto, process, configuration, zero_mass, zero_mass);
		const dis::parameters dis_params{pdf, coupling, process, kinematics, dis_channels, Q2};
		const dis::evaluator dis_evaluator{dis_channels, dis_channels
			| std::views::transform([](const dis::partonic_channel &channel) { return dis::physical_channel{channel}; }) | std::ranges::to<std::flat_set>()};

		const std::vector<sidis::partonic_channel> channels = sidis::channel_factory::create_all_channels(pto, process, configuration, zero_mass, zero_mass);
		const sidis::evaluator evaluator{channels, channels
			| std::views::transform([](const sidis::partonic_channel &channel) { return sidis::physical_channel{channel}; }) | std::ranges::to<std::flat_set>()};

		integration_parameters parameters;
		parameters.method = integration_method::cuba_vegas;
		parameters.cuba.maximum_evaluations = 10'000'000;
		parameters.cuba.maximum_relative_error = 0.0;
		parameters.cuba.vegas.starting_iterations = 1'000'000;
		const integrator_mc integrator{parameters};

		for (const coefficient_implementation_t implementation : {coefficient_implementation_t::hw, coefficient_implementation_t::gc}) {
			const sidis::parameters params{pdf, std::vector{ff}, coupling, process, kinematics, pto, channels, implementation, test_thresholds, Q2, Q2, Q2};

			const auto check = [&]<structure_function_t sf>(const double xi) {
				SCOPED_TRACE(std::format("{}, F{}, xi = {}", implementation, std::to_underlying(sf) + 1, xi));
				const double dis = dis_integrand<sf>(xi, pto, dis_evaluator, dis_params);
				const double scale = std::abs(dis) + std::abs(dis_integrand<sf>(xi, lo, dis_evaluator, dis_params));

				const double integral = integrator.integrate([&](const std::span<const double> point, const std::span<double> components) {
					with_coefficient_functions(pto, implementation, [&](auto &&...coefficient_functions) {
						evaluator.evaluate_observable(xi, point[0], point[1], params, decays, coefficient_functions..., structure_function::observable<sf>{},
							[&](const std::size_t, const std::size_t, const std::size_t, const double value) { components[0] += value; });
					});
				}, {0.0, 0.0}, {1.0, 1.0}, 1).value();

				EXPECT_NEAR(integral, dis, 1e-3 * scale);
			};
			for (const double xi : {0.15, 0.3, 0.5, 0.8, 0.95}) {
				check.template operator()<f1>(xi);
				check.template operator()<f2>(xi);
				check.template operator()<f3>(xi);
			}
		}
	}
}

TEST(momentum_sum_rule, lo) {
	expect_momentum_sum_rule(lo);
}

TEST(momentum_sum_rule, nlo) {
	expect_momentum_sum_rule(nlo);
}
