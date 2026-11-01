#ifndef DM_INTEGRATION_IS_1D_INTEGRATOR_H
#define DM_INTEGRATION_IS_1D_INTEGRATOR_H

#include <concepts>

#include "dimuon/integration/integrator_cquad.hpp"

namespace dm {

/// A one-dimensional integrator should integrate a function of one variable from `lower` to `upper` and return a result object
template <typename integrator_t>
concept is_1d_integrator = requires(const integrator_t integrator, double (*integrand)(double), const double lower, const double upper) {
	{ integrator.integrate(integrand, lower, upper).value } -> std::convertible_to<double>;
};

static_assert(is_1d_integrator<integrator_cquad>, "Type `integrator_cquad` doesn't satisfy concept `is_1d_integrator`.");

}

#endif
