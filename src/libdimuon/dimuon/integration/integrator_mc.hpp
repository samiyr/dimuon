#ifndef DM_INTEGRATION_INTEGRATOR_MC_H
#define DM_INTEGRATION_INTEGRATOR_MC_H

#include <vector>
#include <cassert>
#include <ranges>
#include <utility>
#include <cstddef>
#include <span>
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>

#include <cuba.h>
#include <gsl/gsl_qrng.h>

#include "dimuon/integration/integration_parameters.hpp"
#include "dimuon/integration/integration_result.hpp"

namespace dm {

/// The integrand takes in an input point and a list of components.
/// The integrand is responsible for updating the components with something like `components[component] += value`.
/// The integral is calculated over each component and the total sum of the components.
/// The total is calculated from the totals once all values have been added to avoid chaining the additions.
template <typename integrand_t>
concept is_mc_integrand_function = requires(integrand_t function, const std::span<const double> input, const std::span<double> components) {
	function(input, components);
};

struct integrator_mc {
	/// The purpose of this functor is to enable passing arbitrary callable types to the integration routine.
	/// By default, the Cuba integration library takes a function pointer as the integrand.
	/// Non-capturing lambdas can be implicitly converted to a function pointer, but capturing lambdas
	/// cannot be. To support capturing lambdas, this trick with a functor is used. The idea is to pass a functor
	/// instance as the arbitrary parameter void* to the integration routine. The actual integration parameters are then
	/// stored inside the functor, alongside the actual integrand function and potentially other data as well.
	/// At the integration callsite, instead of passing the actual integrand, the invoke method is passed instead.
	/// Inside the invoke method, the functor instance previously passed as the arbitrary parameter of type void*
	/// is obtained. Using the functor instance the actual integrand function stored in it can then be called.
	template <is_mc_integrand_function integrand_t>
	struct cuba_functor {
		/// Integrand function
		integrand_t function;
		/// An array of lower bounds of the integral, length must be equal to the dimension of the integral
		const std::vector<double> lower;
		/// An array of upper bounds of the integral, length must be equal to the dimension of the integral
		const std::vector<double> upper;
		/// The Jacobian of the integration bounds, i.e. the volume of the integration region
		const double jacobian;

		/// @brief Constructor that also initializes internal variables
		/// @param[in] function_ The integrand function, must satisfy the concept `is_mc_integrand_function`
		/// @param[in] lower_ List of lower integration bounds for each dimension. 
		/// Order must be the same as the order of integration variables.
		/// Length must be equal to the length of \p upper_.
		/// @param[in] upper_ List of upper integration bounds for each dimension. 
		/// Order must be the same as the order of integration variables.
		/// Length must be equal to the length of \p lower_.
		explicit cuba_functor(integrand_t function_, std::vector<double> lower_, std::vector<double> upper_)
		: function{std::move(function_)}, 
		lower(std::move(lower_)), 
		upper(std::move(upper_)), 
		jacobian{cuba_functor::compute_jacobian(lower, upper)}, 
		input(lower.size(), 0.0)
		{ }

		/// @brief Evaluates the integrand at \p x in the unit hypercube, mapped to the integration region, adding its values to \p components
		/// @param[in] x The integration point in the unit hypercube
		/// @param[in,out] components The components, to which the values of the integrand are added
		void evaluate(const double x[], const std::span<double> components) {
			// The integration is always done in a unit hypercube. This loop maps the input
			// to the actual input vector, which may not lie in the unit hypercube.
			for (std::size_t i = 0; i < input.size(); i++) {
				const double range = upper[i] - lower[i];

				input[i] = lower[i] + x[i] * range;
			}

			// Call the actual integrand function with the appropriate parameters
			function(std::span<const double>{input}, components);
		}

		/// @brief This static method is the one that should be passed as a function pointer to the integration routines.
		/// This method will then handle calling the actual integrand and passing any arbitrary `void *` parameters.
		/// The signature of this function is the signature required by Cuba.
		/// @param[in] ndim Dimension of the integral (unused by this function, which uses the dimension of the bounds)
		/// @param[in] x The `ndim`-dimensional integration point, represented as an array
		/// @param[in] ncomp Number of Cuba components: the total and the components of the integrand
		/// @param[out] f Output vector, in which the results of the integration are collected. 
		/// `f[0]` is the total, the sum of the components, and the rest of \p f are the components of the integrand.
		static int invoke(const int */*ndim*/, const double x[], const int *ncomp, double f[], void *userdata) {
			// Obtain the functor passed onto the integration routine as an arbitrary void* parameter
			cuba_functor &functor = *static_cast<cuba_functor *>(userdata);

			// The values of the components of the integrand at this point, which follow the total in the output vector
			const std::span<double> components{f + 1, static_cast<std::size_t>(*ncomp) - 1};
			std::ranges::fill(components, 0.0);
			functor.evaluate(x, components);

			// Fill the output vector
			f[0] = functor.jacobian * std::ranges::fold_left(components, 0.0, std::plus{});
			for (double &component : components) {
				component *= functor.jacobian;
			}

			return 0;
		}

	private:
		/// Underlying input array, which is used to map the unit hypercube into the actual integration region. See `cuba_functor::invoke()`.
		std::vector<double> input;

		/// @brief Compute the Jacobian (i.e. the integration region volume)
		/// @param[in] lower List of lower integration bounds for each dimension. 
		/// Order must be the same as the order of integration variables.
		/// Length must be equal to the length of \p upper.
		/// @param[in] upper List of upper integration bounds for each dimension. 
		/// Order must be the same as the order of integration variables.
		/// Length must be equal to the length of \p lower.
		static double compute_jacobian(const std::vector<double> &lower, const std::vector<double> &upper) noexcept {
			return std::ranges::fold_left(
				std::views::zip_transform([](const double l, const double u) { return u - l; }, lower, upper),
				1.0, std::multiplies{}
			);
		}

	};

	/// Integration result, see `integration_result`
	using result = integration_result;

	/// The method of `integrate`: Cuba's Vegas or Suave, or plain quasi-Monte Carlo
	integration_method method;
	integration_parameters::cuba_parameters cuba;

	explicit integrator_mc(integration_parameters parameters_) noexcept : method{parameters_.method}, cuba{std::move(parameters_.cuba)} {
		disable_cuba_parallelization();
	}

	/// @brief Integrates with the method `method`. Cuba (Vegas or Suave) integrates the total as its first component,
	/// followed by the components of the integrand.
	/// Plain quasi-Monte Carlo takes the mean value of the integrand at `cuba.maximum_evaluations` Sobol points (from GSL)
	/// times the volume of the integration region. The points are those of the first iteration of Cuba's Vegas, and there is no adaptation.
	/// The values are summed as with Cuba, but over all points at once, so the cost of a value does not depend on the number of components.
	/// This suits integrands with many components of which each point has only a few non-zero.
	/// There are no error estimates, since the variance of the values does not estimate the error of quasi-Monte Carlo integration:
	/// the errors and probabilities of the result are NaN.
	/// @param[in] integrand The integrand function, must satisfy the concept `is_mc_integrand_function`
	/// @param[in] lower List of lower integration bounds for each dimension
	/// @param[in] upper List of upper integration bounds for each dimension
	/// @param[in] n_components Number of components of the integrand
	template <is_mc_integrand_function integrand_t>
	result integrate(integrand_t &&integrand, const std::vector<double> &lower, const std::vector<double> &upper, const std::size_t n_components) const {
		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const int dim = static_cast<int>(lower.size());
		// The total, followed by the components
		const std::size_t true_ncomps = n_components + 1;

		std::vector<double> integrals(true_ncomps);
		std::vector<double> errors(true_ncomps);
		std::vector<double> probabilities(true_ncomps);

		cuba_functor<integrand_t> functor{std::forward<integrand_t>(integrand), lower, upper};
		const auto cuba_integrand = &cuba_functor<integrand_t>::invoke;

		long long int neval{0};
		int fail{0};

		switch (method) {
			case integration_method::cuba_vegas:
				llVegas(
					dim, static_cast<int>(true_ncomps), cuba_integrand, &functor, 1, 
					cuba.maximum_relative_error, cuba.maximum_absolute_error,
					cuba.flags(), cuba.seed, cuba.minimum_evaluations, cuba.maximum_evaluations,
					cuba.vegas.starting_iterations, cuba.vegas.continuing_iterations, cuba.vegas.batch_size, cuba.vegas.grid_number,
					nullptr, nullptr, &neval, &fail,
					integrals.data(), errors.data(), probabilities.data()
				);
				break;

			case integration_method::cuba_suave:
			{
				int nregions{0};

				llSuave(
					dim, static_cast<int>(true_ncomps), cuba_integrand, &functor, 1, 
					cuba.maximum_relative_error, cuba.maximum_absolute_error,
					cuba.flags(), cuba.seed, cuba.minimum_evaluations, cuba.maximum_evaluations,
					cuba.suave.new_subdivision_evaluations, cuba.suave.minimum_samples, cuba.suave.flatness,
					nullptr, nullptr, &nregions, &neval, &fail,
					integrals.data(), errors.data(), probabilities.data()
				);

				break;
			}

			case integration_method::plain_qmc:
			{
				const std::unique_ptr<gsl_qrng, decltype(&gsl_qrng_free)> sobol{gsl_qrng_alloc(gsl_qrng_sobol, static_cast<unsigned int>(lower.size())), gsl_qrng_free};
				if (!sobol) {
					throw std::runtime_error("Could not allocate the Sobol generator for quasi-Monte Carlo integration");
				}

				// The integrals start as the sums of the values over all points
				const std::span<double> components = std::span{integrals}.subspan(1);
				std::vector<double> x(lower.size());
				for (long long int point = 0; point < cuba.maximum_evaluations; ++point) {
					gsl_qrng_get(sobol.get(), x.data());
					functor.evaluate(x.data(), components);
				}
				integrals[0] = std::ranges::fold_left(components, 0.0, std::plus{});

				// Each integral is the mean value times the volume of the region
				std::ranges::transform(integrals, integrals.begin(), [&](const double sum) { return functor.jacobian * sum / static_cast<double>(cuba.maximum_evaluations); });
				std::ranges::fill(errors, std::numeric_limits<double>::quiet_NaN());
				std::ranges::fill(probabilities, std::numeric_limits<double>::quiet_NaN());
				neval = cuba.maximum_evaluations;

				break;
			}
		}

		return result {
			integrals, errors, probabilities,
			neval
		};
	}

private:
	void set_cuba_cores(const int cores, const int max) noexcept {
		cubacores(&cores, &max);
	}
	void disable_cuba_parallelization() noexcept {
		if (!cuba.parallelize) {
			set_cuba_cores(0, 0);
		}
	}
};

}

#endif