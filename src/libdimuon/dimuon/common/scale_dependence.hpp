#ifndef DM_COMMON_SCALE_DEPENDENCE_H
#define DM_COMMON_SCALE_DEPENDENCE_H

#include <concepts>
#include <algorithm>
#include <type_traits>

#include "dimuon/common/trf_kinematics.hpp"

namespace dm {

/// Scale dependence concept. Checks that the object can be called with a `trf_kinematics` object and that it returns a double.
template <typename T>
concept is_scale_dependence = requires(T scale_function, const trf_kinematics &kinematics) {
	{ scale_function(kinematics) } -> std::convertible_to<double>;
};

namespace scale_dependence {
	/// A wrapper for a scale dependence function, where the template parameter determines the function type. The function type
	/// must satisfy the 'is_scale_dependence' concept; see is_scale_dependence. The function can be left undefined, in which case
	/// calling the wrapper defaults to returning the Q^2 value.
	template <is_scale_dependence function_t>
	struct function {
		double minimum;

		/// Constructs an empty function wrapper, in which case calling the constructed object returns Q^2
		explicit function(double minimum_ = 0.0) noexcept : minimum{minimum_}, scale_function{[](const trf_kinematics &kinematics) { return kinematics.Q2; }} {}
		/// Constructs a function wrapper, in which case calling the constructed object calls the provided function
		explicit function(function_t scale_function_, double minimum_ = 0.0) noexcept(std::is_nothrow_move_constructible_v<function_t>)
		: minimum{minimum_}, scale_function{std::move(scale_function_)} {}

		/// Evaluates the wrapped function if one exists, otherwise returns Q^2
		double operator()(const trf_kinematics &kinematics) const noexcept(noexcept(scale_function(kinematics))) {
			return std::max(scale_function(kinematics), minimum);
		}

	private:
		function_t scale_function;
	};

	/// Scale dependence function wrapper that returns Q^2
	inline auto Q2 = function([](const trf_kinematics &kinematics) noexcept { return kinematics.Q2; });
	/// Scale dependence function wrapper that returns Q^2, with a non-zero minimum
	[[maybe_unused]] inline auto clamped(const double minimum) noexcept {
		return function([](const trf_kinematics &kinematics) noexcept { return kinematics.Q2; }, minimum);
	}
	/// Scale dependence function wrapper that returns a constant value
	[[maybe_unused]] inline auto constant(const double value) noexcept {
		return function([value](const trf_kinematics &) noexcept { return value; });
	}
	/// Scale dependence function wrapper that returns Q^2 multiplied by a constant value
	[[maybe_unused]] inline auto multiplicative(const double factor) noexcept {
		return function([factor](const trf_kinematics &kinematics) noexcept { return kinematics.Q2 * factor; });
	}
	/// Scale dependence function wrapper that returns Q^2 multiplied by a constant value, in such a way that the multiplied value is clamped
	/// between a minimum and maximum Q^2 value
	[[maybe_unused]] inline auto clamped_multiplicative(const double factor, const double Q2_min, const double Q2_max) noexcept {
		return function([factor, Q2_min, Q2_max](const trf_kinematics &kinematics) noexcept { return std::clamp(kinematics.Q2 * factor, Q2_min, Q2_max); });
	}

	template <is_scale_dependence function_t>
	auto operator*(const double lhs, const function<function_t> &rhs) {
		return function([lhs, rhs](const trf_kinematics &kinematics) noexcept(noexcept(rhs(kinematics))) {
			return lhs * rhs(kinematics);
		}, rhs.minimum);
	}
	
	template <is_scale_dependence function_t>
	auto operator*(const function<function_t> &lhs, const double rhs) {
		return rhs * lhs;
	}
}

}

#endif