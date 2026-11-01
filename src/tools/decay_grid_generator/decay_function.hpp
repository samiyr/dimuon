#ifndef DM_TOOLS_DECAY_GRID_GENERATOR_DECAY_FUNCTION_H
#define DM_TOOLS_DECAY_GRID_GENERATOR_DECAY_FUNCTION_H

#include <concepts>
#include <numbers>
#include <filesystem>

#include <dimuon/common/particle.hpp>

#include <dimuon/utility/math.hpp>

#include <dimuon/decay/decay_parametrization.hpp>

#include <dimuon/interpolation/grid_interpolation.hpp>

inline constexpr double differential_decay_function_integrand(
	const double rho, const double reduced_rho, const double cos, const double a, const double b, const double h0,
	const dm::decay_parametrization &parametrization, const dm::particle &resonance
) {
	const double w = a * rho - b * reduced_rho * cos;

	if (w < 0.0 || w > 1.0 / parametrization.gamma) { return 0.0; }

	const double decay_value = parametrization.N * std::pow(w, parametrization.alpha) * std::pow(1.0 - parametrization.gamma * w, parametrization.beta);
	const double integrand = 2.0 * std::numbers::pi * resonance.lifetime * decay_value * reduced_rho * dm::math::power<2>(h0) / (2.0 * resonance.mass);

	return integrand;
}

#endif