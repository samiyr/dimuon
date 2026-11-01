#include <cmath>
#include <numbers>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/integration/integrator_cquad.hpp>
#include <dimuon/integration/integrator_mc.hpp>
#include <dimuon/integration/integrator_nd.hpp>
#include <dimuon/utility/math.hpp>

#include "support/assertions.hpp"

using namespace dm;

TEST(integrator_cquad, integration) {
	const integrator_cquad integrator{0.0, 1e-12};

	EXPECT_REL_NEAR(integrator.integrate([](const double x) { return math::power<2>(x) * std::exp(-x); }, -1.0, 2.0).value, std::numbers::e - 10.0 / math::power<2>(std::numbers::e), 1e-12);
	EXPECT_REL_NEAR(integrator.integrate([](const double x) { return std::log(x); }, 0.0, 1.0).value, -1.0, 1e-12);
	EXPECT_REL_NEAR(integrator.integrate([](const double x) { return 1.0 / std::sqrt(x); }, 0.0, 1.0).value, 2.0, 1e-12);
}

TEST(integrator_nd, integration) {
	const integrator_nd<integrator_cquad, 3> integrator{0.0, 1e-12, std::size_t{1000}};
	EXPECT_REL_NEAR(integrator.integrate([](const double x, const double y, const double z) { return x * math::power<2>(y) * math::power<3>(z); }, {0.0, 1.0, 0.0}, {1.0, 2.0, 2.0}), 14.0 / 3.0, 1e-10);

	const integrator_nd<integrator_cquad, 2> integrator_2d{0.0, 1e-12, std::size_t{1000}};
	EXPECT_REL_NEAR(integrator_2d.integrate([](const double x, const double y) { return std::exp(x * y); }, {0.0, 0.0}, {1.0, 1.0}), 1.317902151454404, 1e-10);
}

TEST(integrator_mc, vector_integration) {
	const auto integrand = [](const std::span<const double> point, const std::span<double> components) {
		const double x = point[0];
		const double y = point[1];
		const double z = point[2];

		components[0] += math::power<2>(x);
		components[1] += y * z;
	};
	const std::vector<double> lower{0.0, 0.0, 1.0};
	const std::vector<double> upper{2.0, 1.0, 3.0};
	// The integrals of the components over the box
	const double x2_integral = 16.0 / 3.0;
	const double yz_integral = 4.0;

	for (const integration_method method : {integration_method::cuba_vegas, integration_method::cuba_suave, integration_method::plain_qmc}) {
		SCOPED_TRACE(names::name_of(integration_method_names, method));
		integration_parameters parameters;
		parameters.method = method;
		parameters.cuba.maximum_relative_error = 1e-4;
		parameters.cuba.maximum_evaluations = 1'000'000;

		const integrator_mc::result result = integrator_mc{parameters}.integrate(integrand, lower, upper, 2);
		const std::vector<double> components(result.component_values().begin(), result.component_values().end());

		ASSERT_EQ(components.size(), 2);
		EXPECT_REL_NEAR(components[0], x2_integral, 1e-3);
		EXPECT_REL_NEAR(components[1], yz_integral, 1e-3);
		EXPECT_REL_NEAR(result.value(), x2_integral + yz_integral, 1e-3);

		if (method == integration_method::plain_qmc) {
			EXPECT_REL_NEAR(result.value(), components[0] + components[1], 1e-14);
		} else {
			EXPECT_NEAR(result.value(), components[0] + components[1], 5.0 * result.error());
		}
	}
}
