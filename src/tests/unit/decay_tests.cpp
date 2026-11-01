#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include <gtest/gtest.h>

#include <dimuon/common/constants.hpp>
#include <dimuon/common/particle.hpp>
#include <dimuon/common/trf_kinematics.hpp>
#include <dimuon/decay/decay_grid.hpp>
#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/utility/math.hpp>

#include "support/assertions.hpp"

using namespace dm;

namespace {
	/// Polynomial interpolation reproduces a cubic exactly (within floating-point error)
	double polynomial_grid_function(const double zyE) {
		const double l = std::log(zyE);
		return 0.1 + 0.2 * l - 0.03 * math::power<2>(l) + 0.004 * math::power<3>(l);
	}

	/// A non-polynomial function to interpolate
	double transcendental_grid_function(const double zyE) {
		const double l = std::log(zyE);
		return 0.1 + 0.2 * std::sin(l) - 0.03 * std::exp(-math::power<2>(l));
	}

	/// The range of the grids
	const double E_min = 5.0;
	const double E_max = 300.0;

	/// Unevenly spaced knots of the polynomial grids
	const std::vector<double> polynomial_knots{E_min, 8.0, 15.0, 40.0, 100.0, E_max};

	/// \p count knots at the Chebyshev points of log E, like those of the decay grids of the dimuon program. The grids interpolate in log E with a single polynomial,
	/// which converges at the Chebyshev points of log E: 30 of them reproduce `transcendental_grid_function` within about 1e-13. More knots don't help:
	/// with GSL's polynomial interpolation the rounding errors grow above about 40 knots, and with knots at the Chebyshev points of E the interpolation diverges.
	std::vector<double> log_chebyshev_knots(const std::size_t count) {
		std::vector<double> knots = math::chebyshev_space(std::log(E_min), std::log(E_max), count);
		for (double &knot : knots) { knot = std::exp(knot); }
		// Exactly the ends of the range, which exp(log(E)) can miss by a rounding error
		knots.front() = E_min;
		knots.back() = E_max;
		return knots;
	}

	/// Points inside the range of the grids, including its ends
	const std::vector<double> inside{E_min, 6.5, 50.0, 299.0, E_max};

	/// A temporary directory for the grid files of a test, removed with it
	class decay_grid_file : public ::testing::Test {
	protected:
		std::filesystem::path directory;

		void SetUp() override {
			// Unique for each test and process, since the tests may run in parallel
			const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
			directory = std::filesystem::temp_directory_path() / std::format("dimuon-decay-grid-test-{}-{}", getpid(), test_name);
			std::filesystem::create_directories(directory);
		}
		void TearDown() override { std::filesystem::remove_all(directory); }

		/// Writes a grid file of \p grid_function on the logarithms of \p Es, with the version line \p version_line
		std::filesystem::path write_grid(const std::vector<double> &Es, auto &&grid_function, const std::string &version_line = std::format("version: {}", decay_grid::version)) const {
			const std::filesystem::path path = directory / "grid.dat";
			std::ofstream file(path);
			file << version_line << "\nresonance: test\n-----\n";
			for (const double E : Es) { file << std::format("{:.17g}\n", std::log(E)); }
			file << "-----\n";
			for (const double E : Es) { file << std::format("{:.17g}\n", std::forward<decltype(grid_function)>(grid_function)(E)); }
			return path;
		}
	};
}

TEST_F(decay_grid_file, polynomial_interpolation) {
	const decay_grid grid{write_grid(polynomial_knots, polynomial_grid_function)};
	for (const double zyE : inside) {
		EXPECT_REL_NEAR(grid(zyE), polynomial_grid_function(zyE), 1e-12) << "zyE = " << zyE;
	}
}

TEST_F(decay_grid_file, transcendental_interpolation) {
	const decay_grid grid{write_grid(log_chebyshev_knots(30), transcendental_grid_function)};
	for (const double zyE : inside) {
		EXPECT_REL_NEAR(grid(zyE), transcendental_grid_function(zyE), 1e-12) << "zyE = " << zyE;
	}
}

TEST_F(decay_grid_file, out_of_range) {
	const decay_grid grid{write_grid(polynomial_knots, polynomial_grid_function)};
	
	EXPECT_THROW(static_cast<void>(grid(E_min - 0.1)), std::runtime_error);
	EXPECT_THROW(static_cast<void>(grid(E_max + 1.0)), std::runtime_error);
}

TEST_F(decay_grid_file, missing_file) {
	EXPECT_THROW(decay_grid{directory / "missing.dat"}, std::runtime_error);
}

TEST_F(decay_grid_file, invalid_version) {
	EXPECT_THROW(decay_grid{write_grid({5.0, 10.0, 20.0}, polynomial_grid_function, std::format("version: {}", decay_grid::version + 1))}, std::runtime_error);
	EXPECT_THROW(decay_grid{write_grid({5.0, 10.0, 20.0}, polynomial_grid_function, "5.0")}, std::runtime_error);
}

TEST(decay_instance, caching) {
	// Check that caching works and doesn't become stale
	const auto decay_function = [](const double x, const double z, const double Q2, const particle &) { return x + 10.0 * z + 100.0 * Q2; };
	const particle &target = constants::particles::proton;
	const decay_instance decay{constants::particles::D0, target, decay_function, 5.0};
	const auto expect = [&](const double x, const double z, const double Q2) { EXPECT_EQ(decay(x, z, Q2), decay_function(x, z, Q2, target)); };

	const double x1 = 0.1, x2 = 0.3;
	const double z1 = 0.2, z2 = 0.4;
	const double Q2_1 = 3.0, Q2_2 = 5.0;
	expect(x1, z1, Q2_1);
	expect(x2, z1, Q2_1);
	expect(x2, z2, Q2_1);
	expect(x2, z2, Q2_2);
	expect(x1, z1, Q2_1);
}
