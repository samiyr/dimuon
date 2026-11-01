#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/common/structure_function.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/pdf/interface/functional_form_interface.hpp>
#include <dimuon/sidis/channel_factory.hpp>
#include <dimuon/utility/math.hpp>

#include "support/reference_values.hpp"
#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum sidis::partonic_channel::channel_type;
	using enum structure_function_t;

	struct coefficient_case {
		std::string name;
		perturbative_order pto;
		sidis::partonic_channel::channel_type type;
	};

	void PrintTo(const coefficient_case &coefficient, std::ostream *stream) { *stream << coefficient.name; }

	class sidis_reference : public ::testing::TestWithParam<coefficient_case> {};

	struct unit_coupling {
		double operator()(const flavor::type, const flavor::type) const noexcept { return 1.0; }
	};

	const functional_form_interface reference_pdf{[](const flavor::type, const double x, const double) { return math::power<3>(1.0 - x) / (x * x); }};
	const functional_form_interface reference_ff{[](const flavor::type, const double z, const double) { return math::power<3>(z) * std::pow(1.0 - z, 1.5); }};

	structure_function_point make_point(const double x, const double z) { return {x, z, 4.0, 4.0, 4.0, 4.0}; }

	std::string order_name(const perturbative_order pto) {
		switch (pto) {
			case perturbative_order::lo: return "lo";
			case perturbative_order::nlo: return "nlo";
			case perturbative_order::nnlo: return "nnlo";
		}
		std::unreachable();
	}

	std::string reference_file(const structure_function_t sf, const coefficient_case &coefficient) {
		return std::format("sidis/F{}{}.tsv", std::to_underlying(sf) + 1, coefficient.name);
	}

	sidis::partonic_channel make_channel(const coefficient_case &coefficient) {
		const scattering_process process{scattering_process::process_type::neutrino, constants::particles::proton, constants::particles::neutrino};
		const channel_configuration zero_mass{mass_scheme_t::zero_mass};
		const std::vector<sidis::partonic_channel> all = sidis::channel_factory::create_all_channels(
			coefficient.pto, process, parton_configuration{std::vector{flavor::s, flavor::c}, std::array<double, constants::total_flavors>{}}, zero_mass, zero_mass
		);
		const auto channel = std::ranges::find(all, coefficient.type, &sidis::partonic_channel::type);
		if (channel == all.end()) { throw std::runtime_error(std::format("no channel of the coefficient function {}", coefficient.name)); }
		return *channel;
	}

	double compute(const structure_function_t sf, const sidis::partonic_channel &channel, const coefficient_implementation_t implementation, const structure_function_point &point) {
		integration_parameters parameters;
		parameters.method = integration_method::cuba_vegas;
		parameters.cuba.maximum_evaluations = 10'000'000;
		parameters.cuba.maximum_relative_error = 0.0;
		parameters.cuba.vegas.starting_iterations = 1'000'000;
		const auto compute_sf = [&]<structure_function_t sf_>() {
			return fixed_z_structure_function<sf_>(std::vector{channel}, channel.pto, implementation, reference_pdf, reference_ff, unit_coupling{}, point, parameters);
		};
		switch (sf) {
			case f1: return point.x * point.z * compute_sf.template operator()<f1>();
			case f2: return point.z * compute_sf.template operator()<f2>();
			case f3: return point.x * point.z * compute_sf.template operator()<f3>();
		}
		std::unreachable();
	}
}

TEST_P(sidis_reference, coefficient_function) {
	const coefficient_case &coefficient = GetParam();
	const sidis::partonic_channel channel = make_channel(coefficient);
	const double tolerance = coefficient.pto == perturbative_order::nnlo ? 5e-3 : 5e-5;
	std::vector<coefficient_implementation_t> implementations{coefficient_implementation_t::gc};
	if (coefficient.pto != perturbative_order::nnlo) { implementations.push_back(coefficient_implementation_t::hw); }

	const std::string f1_file = reference_file(f1, coefficient);
	const std::optional<std::vector<reference_point>> f1_reference = read_reference_points(f1_file);
	if (!f1_reference) { GTEST_SKIP() << missing_reference_message(f1_file); }
	const auto f1_value = [&](const std::vector<double> &inputs) {
		const auto row = std::ranges::find(*f1_reference, inputs, &reference_point::inputs);
		return row == f1_reference->end() ? std::nullopt : std::optional{row->output};
	};

	for (const structure_function_t sf : {f1, f2, f3}) {
		const std::string file = reference_file(sf, coefficient);
		const std::optional<std::vector<reference_point>> reference = read_reference_points(file);

		// The files of the coefficient functions whose F3 is zero have only F1 and F2
		if (!reference) {
			EXPECT_EQ(sf, f3) << missing_reference_message(file);
			continue;
		}
		for (const auto &[inputs, expected] : *reference) {
			ASSERT_EQ(inputs.size(), 2) << file << ": a row is not \"x z value\"";

			const double x = inputs[0];
			const double z = inputs[1];
			double scale = std::abs(expected);

			if (sf == f2) {
				const std::optional<double> f1_expected = f1_value(inputs);
				ASSERT_TRUE(f1_expected) << f1_file << " has no row at x = " << x << ", z = " << z;
				scale += 2.0 * std::abs(*f1_expected);
			}
			for (const coefficient_implementation_t implementation : implementations) {
				EXPECT_NEAR(compute(sf, channel, implementation, make_point(x, z)), expected, tolerance * scale)
					<< names::name_of(coefficient_implementation_names, implementation) << ", " << file << ", x = " << x << ", z = " << z;
			}
		}
	}
}

INSTANTIATE_TEST_SUITE_P(
	coefficients, sidis_reference,
	::testing::Values(
		coefficient_case{"C0q2qM", perturbative_order::lo, beta_alpha},
		coefficient_case{"C1q2qM", perturbative_order::nlo, beta_alpha},
		coefficient_case{"C1q2gM", perturbative_order::nlo, g_alpha},
		coefficient_case{"C1g2qM", perturbative_order::nlo, beta_g},
		coefficient_case{"C2q2qM", perturbative_order::nnlo, beta_alpha},
		coefficient_case{"C2q2qMFcon1", perturbative_order::nnlo, alpha_alpha_Fcon},
		coefficient_case{"C2q2qMFcon2", perturbative_order::nnlo, qq_Fcon},
		coefficient_case{"C2q2qbM", perturbative_order::nnlo, alpha_betabar},
		coefficient_case{"C2q2qbMFcon", perturbative_order::nnlo, alphabar_alpha_Fcon},
		coefficient_case{"C2q2qpM1", perturbative_order::nnlo, qprime_alpha},
		coefficient_case{"C2q2qpM2", perturbative_order::nnlo, betaprime_q},
		coefficient_case{"C2q2gM", perturbative_order::nnlo, g_alpha},
		coefficient_case{"C2g2qM", perturbative_order::nnlo, beta_g},
		coefficient_case{"C2g2gM", perturbative_order::nnlo, gg}
	),
	[](const auto &parameter) { return std::format("{}_{}", order_name(parameter.param.pto), parameter.param.name); }
);
