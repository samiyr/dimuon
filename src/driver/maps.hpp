#ifndef DM_EXE_MAPS_H
#define DM_EXE_MAPS_H

#include <cstddef>
#include <ranges>
#include <unordered_map>
#include <string>
#include <utility>

#include <dimuon/common/particle.hpp>
#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/mass_scheme.hpp>
#include <dimuon/common/coefficient_implementation.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/utility/names.hpp>

#include "run_mode.hpp"
#include "output_mode.hpp"

namespace maps {
	/// Returns a map from the bidirectional map (see `dm::names`)
	template <typename value_t, std::size_t size>
	std::unordered_map<std::string, value_t> by_name(const dm::names::table<value_t, size> &names) {
		return names
			| std::views::transform([](const dm::names::entry<value_t> &entry) { return std::pair{std::string{entry.second}, entry.first}; })
			| std::ranges::to<std::unordered_map<std::string, value_t>>();
	}

	static const std::unordered_map<std::string, run_mode> run_modes = {
		{"direct", run_mode::direct},
		{"kernel", run_mode::kernel}
	};

	static const std::unordered_map<std::string, output_mode> output_modes = {
		{"stdout", output_mode::standard_output},
		{"file", output_mode::file}
	};

	static const std::unordered_map<std::string, dm::perturbative_order> perturbative_orders = {
		{"lo", dm::perturbative_order::lo},
		{"nlo", dm::perturbative_order::nlo},
		{"nnlo", dm::perturbative_order::nnlo}
	};

	/// The scale dependences (--scale), by whether the squared charm mass is added to Q^2: mu^2 = Q^2 (Q2) or mu^2 = Q^2 + mc^2 (Q2+mc2)
	static const std::unordered_map<std::string, bool> scale_dependences = {
		{"Q2", false},
		{"Q2+mc2", true}
	};

	static const std::unordered_map<std::string, dm::mass_scheme_t> mass_schemes = by_name(dm::mass_scheme_names);

	static const std::unordered_map<std::string, dm::coefficient_implementation_t> coefficient_implementations = by_name(dm::coefficient_implementation_names);

	static const std::unordered_map<std::string, dm::integration_method> integration_methods = by_name(dm::integration_method_names);

	static const std::unordered_map<std::string, dm::particle> resonance_particles = {
		{"D0", dm::constants::particles::D0},
		{"Dp", dm::constants::particles::Dp},
		{"Ds", dm::constants::particles::Ds},
		{"LambdaC", dm::constants::particles::LambdaC}
	};
}

#endif