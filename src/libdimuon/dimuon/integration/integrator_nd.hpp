#ifndef DM_INTEGRATION_INTEGRATOR_ND_H
#define DM_INTEGRATION_INTEGRATOR_ND_H

#include <cstddef>
#include <array>
#include <concepts>
#include <utility>

#include "dimuon/integration/is_1d_integrator.hpp"

namespace dm {

/// @brief Integrates a function of `dim` variables by nesting one-dimensional integrals: the integrand of the integral over each variable
/// is the integral over the variables after it. Each variable has its own one-dimensional integrator, since its integrals are computed
/// while the integral over the previous variable is in progress. The one-dimensional integrators may keep state between their integrals
/// (such as the workspace of `integrator_cquad`), so an instance must not be used by several threads at once, or from its own integrand.
/// @tparam integrator_1d_t The type of the one-dimensional integrators, see `is_1d_integrator`
/// @tparam dim The number of variables
template <is_1d_integrator integrator_1d_t, std::size_t dim>
struct integrator_nd {
	/// Constructs the one-dimensional integrator of each variable from \p args
	template <typename... args_t>
	requires std::constructible_from<integrator_1d_t, const args_t &...>
	explicit integrator_nd(const args_t &...args) : integrators{construct_integrators(std::make_index_sequence<dim>{}, args...)} {}

	template <typename integrand_t>
	double integrate(
		integrand_t &&integrand,
		const std::array<double, dim> &lower, const std::array<double, dim> &upper
	) const {
		return recurse(std::forward<integrand_t>(integrand), lower, upper);
	}

private:
	/// The one-dimensional integrator of each variable
	std::array<integrator_1d_t, dim> integrators;

	/// Returns a one-dimensional integrator constructed from \p args for each variable
	template <std::size_t... variable, typename... args_t>
	static std::array<integrator_1d_t, dim> construct_integrators(std::index_sequence<variable...>, const args_t &...args) {
		return {(static_cast<void>(variable), integrator_1d_t(args...))...};
	}

	/// This function does the actual integration by recursion. The constant template parameter `n`,
	/// which defaults to `0`, denotes the depth of the recursion. `0` depth means no recursion yet.
	/// While `n < dim - 1`, this function will recursively call itself as the integrand, while
	/// every time incrementing `n`. The parameter pack `accumulated_input` collects the integration variables.
	/// During the first call (`n = 0`), the `accumulated_input` pack is empty. The integration routine will then
	/// provide values for the integration variables one level of recursion at a time, which then starts
	/// to fill the `accumulated_input` pack. The recursion is broken when `n == dim - 1` by actually calling the
	/// provided `integrand` function with the fully-collected `accumulated_input` pack.
	///
	/// For example, a 3D integral would look like this 
	/// (with `lower` and `upper` implicit):
	/// ```
	/// recurse<0>(f(., ., .))
	///	--> integrate(
	///			recurse<1>(f(x, ., .)), 
	///		x)
	///	--> integrate(
	///			integrate(
	///				recurse<2>(f(x, y, .)), 
	///			y),
	///		x)
	///	--> integrate(
	///			integrate(
	///				integrate(
	///					f(x, y, z)
	///				z)
	///			y)
	///		x)
	/// ```
	template <std::size_t n = 0, typename integrand_t, std::floating_point... input_t>
	double recurse(
		integrand_t &&integrand,
		const std::array<double, dim> &lower, const std::array<double, dim> &upper,
		input_t... accumulated_input
	) const {
		// Integrate over dimension `n` with its own one-dimensional integrator
		return integrators[n].integrate([&](const double current_input) {
			if constexpr (n == dim - 1) {
				// The recursion is broken once the final integration dimension is reached.
				// At this point, the actual integrand function can be called with the collected
				// input, including the variable of the final dimension.
				return std::forward<integrand_t>(integrand)(accumulated_input..., current_input);
			} else {
				// Continue the recursion by appending the variable of the current dimension
				// to the input parameter pack and incrementing the recursion depth.
				return recurse<n + 1>(std::forward<integrand_t>(integrand), lower, upper, accumulated_input..., current_input);
			}
		}, lower[n], upper[n]).value;
		// Note that the order of appending the `current_input` to `accumulated_input` matters, since
		// the integration interval for each dimension is given by `lower[n]` and `upper[n]`.
	}

};

}

#endif