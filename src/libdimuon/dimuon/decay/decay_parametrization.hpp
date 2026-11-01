#ifndef DM_DECAY_PARAMETRIZATION_H
#define DM_DECAY_PARAMETRIZATION_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <format>

namespace dm {

/// Parameters for the decay function functional form
struct decay_parametrization {
	/// Returns the default fit parametrization
	constexpr static decay_parametrization default_fit() noexcept {
		return decay_parametrization{3.1561790889845995, 1.0087668366926346, 1.7993246621566206, 2.100989529372037};
	}

	/// Loads a list of parametrizations from a CSV file, where each line is of the form
	/// `N,alpha,beta,gamma`
	static std::vector<decay_parametrization> fit_set(const std::string &filename) {
		std::ifstream file(filename);
		if (!file) {
			throw std::runtime_error(std::format("Fit set file '{}' could not be opened.", filename));
		}

		std::vector<std::vector<double>> parameter_set;

		std::size_t current_index = 0;

		std::string line;
		while (std::getline(file, line)) {
			current_index++;

			if (current_index == 1) { continue; }

			std::vector<double> parameters;

			std::istringstream stream(line);
			std::string value;

			while (std::getline(stream, value, ',')) {
				parameters.push_back(std::stod(value));
			}

			parameter_set.push_back(parameters);
		}

		std::vector<decay_parametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}

	double N{};
	double alpha{};
	double beta{};
	double gamma{};
};

}

#endif