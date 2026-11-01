#ifndef DM_DIS_COEFFICIENTS_COEFFICIENT_FUNCTION_H
#define DM_DIS_COEFFICIENTS_COEFFICIENT_FUNCTION_H

#include <concepts>

#include "dimuon/dis/coefficients/evaluation_parameters.hpp"
#include "dimuon/dis/vectorization.hpp"

namespace dm {

namespace dis::coefficients {
	/// A concept of a coefficient function, such as `f1::nlo::quark`: the integrand of one structure function at one perturbative order,
	/// in the quark or gluon channels (see `dis::evaluator::evaluate_observable`).
	/// `function(p)` gives its value with the evaluation parameters `p` (see `dis::evaluation_parameters`), with one lane per PDF value
	/// (see `dis::vectorization::distribution_product`).
	template <typename function_t>
	concept coefficient_function = requires(function_t function, const dis::evaluation_parameters &p) {
		{ function(p) } -> std::convertible_to<dis::vectorization::distribution_product>;
	};
}

}

#endif
