#ifndef DM_INTEGRATION_IS_1D_INTEGRAND_FUNCTION_H
#define DM_INTEGRATION_IS_1D_INTEGRAND_FUNCTION_H

#include <concepts>

namespace dm {

/// A 1D integrand should take in the scalar input point and return a scalar
template <typename integrand_t>
concept is_1d_integrand_function = requires(integrand_t function, const double x) {
	{ function(x) } -> std::convertible_to<double>;
};

}

#endif
