#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gsl/gsl_interp.h>
#include <gtest/gtest.h>

#include <dimuon/common/flavor.hpp>
#include <dimuon/interpolation/grid_interpolation.hpp>
#include <dimuon/interpolation/linear_interpolation_basis.hpp>
#include <dimuon/pdf/interface/vector_kernel_interface.hpp>
#include <dimuon/utility/math.hpp>

#include "support/assertions.hpp"

using namespace dm;

namespace {
	/// Unevenly spaced knots
	const std::vector<double> knots{0.01, 0.05, 0.1, 0.3, 0.32, 0.7, 1.0};

	/// Points between the knots
	const std::vector<double> inside{0.01, 0.02, 0.05, 0.07, 0.3, 0.31, 0.5, 0.99, 1.0};
	/// Points outside the knots
	const std::vector<double> outside{0.0, 0.005, 1.0 + 1e-12, 1.5};

	/// The values of all kernels at \p x, indexed by kernel
	std::vector<double> kernel_values(const linear_interpolation_basis &basis, const double x) {
		std::vector<double> values(basis.size(), 0.0);
		for (const linear_interpolation_basis::kernel &kernel : basis.nonzero_at(x)) { values[kernel.index] = kernel(x); }
		return values;
	}

	/// The sum of f(knot k) times kernel k at \p x
	double interpolate(const linear_interpolation_basis &basis, const auto &f, const double x) {
		const std::vector<double> values = kernel_values(basis, x);
		const auto function_values = knots | std::views::transform(f);
		return std::inner_product(values.begin(), values.end(), function_values.begin(), 0.0);
	}

	double cubic(const double x) { return 1.0 - 2.0 * x + 3.0 * math::power<2>(x) + 4.0 * math::power<3>(x); }

	/// The range of the interpolation of `cubic`
	const double cubic_min = -1.0;
	const double cubic_max = 2.0;

	/// Polynomial interpolation of `cubic`
	grid_interpolation make_cubic_interpolation() {
		const std::vector<double> points{math::chebyshev_space(cubic_min, cubic_max, 10)};
		return grid_interpolation{points, points | std::views::transform(cubic) | std::ranges::to<std::vector>(), gsl_interp_polynomial};
	}
}

TEST(linear_interpolation_basis, size) {
	EXPECT_EQ(linear_interpolation_basis{knots}.size(), knots.size());
}

TEST(linear_interpolation_basis, knots) {
	// phi_i(x_j) = delta_{ij}
	const linear_interpolation_basis basis{knots};
	for (std::size_t j = 0; j < knots.size(); ++j) {
		const std::vector<double> values = kernel_values(basis, knots[j]);
		for (std::size_t i = 0; i < knots.size(); ++i) { EXPECT_EQ(values[i], i == j ? 1.0 : 0.0) << "kernel " << i << " at knot " << j; }
	}
}

TEST(linear_interpolation_basis, partition_of_unity) {
	const linear_interpolation_basis basis{knots};
	for (const double x : inside) { EXPECT_NEAR(std::ranges::fold_left(kernel_values(basis, x), 0.0, std::plus{}), 1.0, 1e-15) << "x = " << x; }
	for (const double x : outside) { EXPECT_TRUE(std::ranges::empty(basis.nonzero_at(x))) << "x = " << x; }
}

TEST(linear_interpolation_basis, locality) {
	// At most two kernels are non-zero, and none is negative
	const linear_interpolation_basis basis{knots};
	for (const double x : inside) {
		EXPECT_LE(std::ranges::distance(basis.nonzero_at(x)), 2) << "x = " << x;
		for (const double value : kernel_values(basis, x)) { EXPECT_GE(value, 0.0) << "x = " << x; }
	}
}

TEST(linear_interpolation_basis, linear_functions) {
	// Linear functions are reproduced exactly
	const linear_interpolation_basis basis{knots};
	const auto f = [](const double x) { return 1.0 - 2.0 * x; };
	for (const double x : inside) { EXPECT_NEAR(interpolate(basis, f, x), f(x), 1e-14) << "x = " << x; }
}

TEST(linear_interpolation_basis, invalid_knots) {
	EXPECT_THROW(linear_interpolation_basis{std::vector{0.5}}, std::invalid_argument);
	EXPECT_THROW((linear_interpolation_basis{std::vector{0.1, 0.3, 0.2}}), std::invalid_argument);
	EXPECT_THROW((linear_interpolation_basis{std::vector{0.1, 0.2, 0.2, 0.3}}), std::invalid_argument);
}

TEST(grid_interpolation, polynomial) {
	const grid_interpolation interpolation = make_cubic_interpolation();
	for (const double x : {cubic_min, -0.7, 0.0, 0.3, 1.9, cubic_max}) { EXPECT_REL_NEAR(interpolation(x), cubic(x), 1e-12) << "x = " << x; }
}

TEST(grid_interpolation, range) {
	const grid_interpolation interpolation = make_cubic_interpolation();
	EXPECT_EQ(interpolation.range(), (std::pair{cubic_min, cubic_max}));

	const double x = 0.3;
	ASSERT_TRUE(interpolation.try_interpolate(x).has_value());
	EXPECT_REL_NEAR(*interpolation.try_interpolate(x), cubic(x), 1e-12);
	EXPECT_FALSE(interpolation.try_interpolate(cubic_min - 0.01).has_value());
	EXPECT_FALSE(interpolation.try_interpolate(cubic_max + 0.01).has_value());
}

TEST(grid_interpolation, clamping) {
	const grid_interpolation interpolation = make_cubic_interpolation();
	EXPECT_REL_NEAR(interpolation.interpolate_clamped(-5.0), cubic(cubic_min), 1e-12);
	EXPECT_REL_NEAR(interpolation.interpolate_clamped(5.0), cubic(cubic_max), 1e-12);
}
