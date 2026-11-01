#ifndef DM_SIDIS_COEFFICIENTS_COEFFICIENT_FUNCTION_H
#define DM_SIDIS_COEFFICIENTS_COEFFICIENT_FUNCTION_H

#include <concepts>

#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"
#include "dimuon/sidis/vectorization.hpp"

namespace dm {

namespace sidis::coefficients {
	/// A concept of a coefficient function, such as `f1_gc::nlo::quark_to_quark`: the integrand of one structure function at one perturbative order,
	/// in the quark-to-quark, quark-to-gluon, gluon-to-quark, or gluon-to-gluon channels (see `sidis::evaluator::evaluate_observable`).
	/// `function(p)` gives its value with the evaluation parameters `p` (see `sidis::evaluation_parameters`), with one lane per product of a PDF value
	/// and an FF value (see `sidis::vectorization::distribution_product`).
	template <typename function_t>
	concept coefficient_function = requires(function_t function, const sidis::evaluation_parameters &p) {
		{ function(p) } -> std::convertible_to<sidis::vectorization::distribution_product>;
	};
}

}

#endif
