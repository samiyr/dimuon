#ifndef DM_INTEGRATION_RESULT_H
#define DM_INTEGRATION_RESULT_H

#include <iostream>
#include <ranges>
#include <utility>
#include <vector>

namespace dm {

/// Integration result
struct integration_result {
	explicit integration_result(
		std::vector<double> values_, 
		std::vector<double> errors_, 
		std::vector<double> probabilities_,
		long long int evaluations_ = -1
	) noexcept
	: values{std::move(values_)},
	errors{std::move(errors_)},
	probabilities{std::move(probabilities_)},
	n_evaluations{evaluations_} {}

	/// Total value of the integral
	double value() const noexcept {
		return values[0];
	}
	/// Error estimate (if available, otherwise `NaN`)
	double error() const noexcept {
		return errors[0];
	}
	/// Estimate of the chi-squared probability (if available, otherwise `NaN`)
	double probability() const noexcept {
		return probabilities[0];
	}
	/// Number of times the integrand was evaluated
	long long int evaluations() const noexcept {
		return n_evaluations;
	}

	/// List of integral component values
	auto component_values() const noexcept {
		return std::ranges::subrange{values.cbegin() + 1, values.cend()};
	}
	/// List of integral error estimates for the components (if available, otherwise `NaN`)
	auto component_errors() const noexcept {
		return std::ranges::subrange{errors.cbegin() + 1, errors.cend()};
	}
	/// List of chi-squared probability estimates for the components (if available, otherwise `NaN`)
	auto component_probabilities() const noexcept {
		return std::ranges::subrange{probabilities.cbegin() + 1, probabilities.cend()};
	}

	friend std::ostream &operator <<(std::ostream &os, const integration_result &r) {
		return os << r.value() << " +- " << r.error() << " (" << r.probability() << ")";
	}

private:
	std::vector<double> values;
	std::vector<double> errors;
	std::vector<double> probabilities;
	long long int n_evaluations;
};

}

#endif
