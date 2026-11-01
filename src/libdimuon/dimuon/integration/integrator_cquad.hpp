#ifndef DM_INTEGRATION_INTEGRATOR_CQUAD_H
#define DM_INTEGRATION_INTEGRATOR_CQUAD_H

#include <iostream>
#include <memory>

#include <gsl/gsl_integration.h>

#include "dimuon/integration/integration_parameters.hpp"
#include "dimuon/integration/is_1d_integrand_function.hpp"

namespace dm {

/// A simple one-dimensional integrator using CQUAD (doubly-adaptive Clenshaw-Curtis quadrature) of the GSL quadrature library
struct integrator_cquad {
	/// See `integrator_mc::cuba_functor`
	template <is_1d_integrand_function integrand_t>
	struct gsl_1d_functor {
		integrand_t function;

		static double invoke(double x, void *p) {
			gsl_1d_functor *functor = static_cast<gsl_1d_functor *>(p);
			return functor->function(x);
		}
	};

	/// Integration result
	struct result {
		/// Integral value
		double value;
		/// Error estimate
		double error;
		/// Number of times the integrand was evaluated
		std::size_t nevals;
		/// GSL return code
		int return_code;

		friend std::ostream &operator <<(std::ostream &os, const result &r) {
			return os << r.value << " +- " << r.error;
		}
	};

	/// Absolute error target
	double epsabs;
	/// Relative error target
	double epsrel;
	/// Number of intervals of the CQUAD workspace (see `gsl_integration_cquad_workspace`)
	const std::size_t interval_count;

	explicit integrator_cquad(double epsabs_ = 0.0, double epsrel_ = 1e-3, std::size_t interval_count_ = 1'000) noexcept
	: epsabs{epsabs_}, epsrel{epsrel_}, interval_count{interval_count_}, workspace{gsl_integration_cquad_workspace_alloc(interval_count_), gsl_workspace_deleter{}} { }

	explicit integrator_cquad(const integration_parameters &parameters) noexcept
	: integrator_cquad(parameters.cquad.max_absolute_error, parameters.cquad.max_relative_error, parameters.cquad.interval_count) {}

	/// Calculates the integral of \p integrand from \p lower to \p upper
	template <is_1d_integrand_function integrand_t>
	result integrate(integrand_t &&integrand, const double lower, const double upper) const {
		double integral = 0.0;
		double error = 0.0;
		std::size_t nevals = 0;

		gsl_1d_functor<integrand_t> functor { std::forward<integrand_t>(integrand) };
		gsl_function function;

		function.function = &gsl_1d_functor<integrand_t>::invoke;
		function.params = &functor;

		const int return_code = gsl_integration_cquad(&function, lower, upper, epsabs, epsrel, workspace.get(), &integral, &error, &nevals);

		return result { integral, error, nevals, return_code };
	}

private:
	struct gsl_workspace_deleter {
		void operator()(gsl_integration_cquad_workspace *ws) noexcept {
			gsl_integration_cquad_workspace_free(ws);
		}
	};
	std::unique_ptr<gsl_integration_cquad_workspace, gsl_workspace_deleter> workspace;
};

}

#endif