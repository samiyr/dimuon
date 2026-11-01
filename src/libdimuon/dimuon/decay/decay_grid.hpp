#ifndef DM_DECAY_GRID_H
#define DM_DECAY_GRID_H

#include <cmath>
#include <string>
#include <format>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gsl/gsl_interp.h>

#include "dimuon/common/particle.hpp"
#include "dimuon/decay/decay_parametrization.hpp"
#include "dimuon/interpolation/grid_interpolation.hpp"

namespace dm {

/// Callable representation of a decay grid. 
/// The decay grids are implemented with polynomial interpolation, with Chebyshev points in log-space.
struct decay_grid {
	/// The grid file format version.
	/// When a grid file is loaded, the version in the grid files is checked against this version.
	/// If the versions don't match, an exception is thrown.
	static constexpr int version = 1;

	/// Loads a decay grid from file
	explicit decay_grid(const std::filesystem::path &grid_path) : interpolation_function{read_grid_file(grid_path), gsl_interp_polynomial} {}

	/// Convenience call operator to the decay grid interpolation
	double operator()(
		const double x, 
		const double z, 
		const double Q2,
		const particle &target
	) const {
		const double zyE = z * Q2 / (2.0 * x * target.mass);
		return operator()(zyE);
	}
	/// @brief Decay grid interpolation
	/// @return The interpolated decay function value at \p zyE
	/// @throws `std::runtime_error` if \p zyE is outside the grid
	double operator()(const double zyE) const {
		const std::optional<double> value = interpolation_function.try_interpolate(std::log(zyE));
		if (!value) {
			const auto [log_zyE_min, log_zyE_max] = interpolation_function.range();
			throw std::runtime_error(std::format(
				"Decay grid: argument zyE = {} out of range [{}, {}]", zyE, std::exp(log_zyE_min), std::exp(log_zyE_max)
			));
		}
		return *value;
	}

	/// Constructs the filename of the requested decay grid
	static std::string filename(const double E_min, const decay_parametrization &parametrization, const particle &resonance, const particle &target, const particle &lepton) {
		return std::format(
			"dg_E_{:.4f}_r_{:.4f}_{:.4f}_t_{:.4f}_l_{:.4f}_p_{:.6f}_{:.6f}_{:.6f}_{:.6f}.dat",
			E_min, resonance.mass, resonance.lifetime, target.mass, lepton.mass, parametrization.N, parametrization.alpha, parametrization.beta, parametrization.gamma
		);
	}

private:
	grid_interpolation interpolation_function;

	/// @brief Returns the grid points and the grid values of the grid file \p grid_path.
	/// The file has three parts, separated by lines "-----":
	/// 1. metadata, with one "key: value" per line; the first line is the version
	/// 2. grid points
	/// 3. grid values. 
	/// @throws std::runtime_error if the file cannot be opened, or if the versions don't match (see `version`)
	static std::pair<std::vector<double>, std::vector<double>> read_grid_file(const std::filesystem::path &grid_path) {
		std::ifstream grid_file(grid_path);
		if (!grid_file) {
			throw std::runtime_error(std::format("Grid file '{}' could not be opened.", grid_path.generic_string()));
		}

		const std::string separator = "-----";
		const std::string version_line = std::format("version: {}", version);

		std::string line;
		std::getline(grid_file, line);
		if (line != version_line) {
			throw std::runtime_error(std::format(
				"In grid file '{}': expected '{}', got '{}' instead",
				grid_path.generic_string(), version_line, line
			));
		}

		while (std::getline(grid_file, line) && line != separator) {}

		std::vector<double> grid_points;
		while (std::getline(grid_file, line) && line != separator) {
			grid_points.push_back(std::stod(line));
		}

		std::vector<double> grid_values;
		while (std::getline(grid_file, line)) {
			grid_values.push_back(std::stod(line));
		}

		return std::make_pair(grid_points, grid_values);
	}
};

}

#endif