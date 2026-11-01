#ifndef DM_UTILITY_MATH_H
#define DM_UTILITY_MATH_H

#include <cstddef>
#include <vector>
#include <cmath>
#include <utility>
#include <numbers>
#include <cfloat> // IWYU pragma: keep (needed by the #include'd polylogarithm sources)
#include <complex> // IWYU pragma: keep (needed by the #include'd polylogarithm sources)
#include <limits> // IWYU pragma: keep (needed by the #include'd polylogarithm sources)

namespace dm {

namespace math {
	/// The external polylogarithm library (see external/polylogarithm/) contains both header and source files. 
	/// The source files are included here as if they were header files, to keep this library header-only.
	/// The unnamed namespace gives their functions internal linkage, which allows several translation units of a program to include this header. 
	/// The standard headers that the files include are included above, because they must stay outside of the namespaces.
	namespace {
		#include <polylogarithm/Li2.cpp>
		#include <polylogarithm/Li3.cpp>
	}

	/// Computes the integer power of an arbitrary type using exponentiation-by-squaring.
	/// The power gets expanded to repeated squaring by multiplication at compile time. For example,
	/// `power<2>(x)` gets expanded to `x * x`, `power<3>(x)` to `power<2>(x) * x`, `power<4>(x)` to `power<2>(power<2>(x))` and so forth.
	template <int n, typename T> requires requires(T x) { { T{1} } -> std::convertible_to<T>; { x * x } -> std::convertible_to<T>; { x / x } -> std::convertible_to<T>; }
	constexpr T power(T x) noexcept {
		if constexpr (n < 0) { return T{1} / power<-n>(x); }
		else if constexpr (n == 0) { return T{1}; }
		else if constexpr (n == 1) { return x; }
		else if constexpr (n == 2) { return x * x; }
		else return power<2>(power<n / 2>(x)) * power<n % 2>(x);
	}

	/// Computes the polylogarithm for `n == 2` and `n == 3` and `x < 1.0`. If `x >= 1.0`, returns `0.0`.
	template <std::size_t n> requires(n == 2 || n == 3)
	inline double polylogT(const double x) noexcept {
		if (x >= 1.0) { return 0.0; }

		if constexpr (n == 2) {
			return polylogarithm::Li2(x);
		} else if constexpr (n == 3) { 
			return polylogarithm::Li3(x);
		}
		std::unreachable();
	}

	/// Computes the inverse tangent integral `Ti_2(x) = int_0^x arctan(t) / t dt` using a rational approximation
	inline double invtanint(const double input) noexcept {
		// Default values for when `0 <= input <= 1`
		double x = input;
		double multiplier = 1.0;
		double remainder = 0.0;

		// Range reduction of the argument to [0, 1].
		if (input < -1.0) {
			x = -1.0 / input;
			multiplier = -1.0;
			remainder = -0.5 * std::numbers::pi * std::log(-input);
		} else if (-1.0 <= input && input < 0.0) {
			x = -input;
			multiplier = -1.0;
			remainder = 0.0;
		} else if (x > 1.0) {
			x = 1.0 / input;
			multiplier = 1.0;
			remainder = 0.5 * std::numbers::pi * std::log(input);
		}

		constexpr static double P[] = {
			1.000000000000000E+00,
			1.359606750543622E+00,
			2.294603108673590E+00,
			1.965431829017859E+00,
			1.479377395374053E+00,
			7.303876012387028E-01,
			2.481653896797473E-01,
			4.618705978228984E-02,
			6.316042475722170E-04
		};
		constexpr static double Q[] = {
			1.000000000000000E+00,
			1.359606750543613E+00,
			2.405714219785837E+00,
			2.116499245692448E+00,
			1.706678976619235E+00,
			9.111698955486797E-01,
			3.619761589327320E-01,
			9.051409249480067E-02,
			9.340084117634812E-03
		};

		const double x2 = math::power<2>(x);
		const double x4 = math::power<4>(x);
		const double x6 = math::power<6>(x);

		const double p = (P[0] + x * P[1]) + x2 * (P[2] + x * P[3]) + x4 * (P[4] + x * P[5]) + x6 * (P[6] + x * P[7] + x2 * P[8]);
		const double q = (Q[0] + x * Q[1]) + x2 * (Q[2] + x * Q[3]) + x4 * (Q[4] + x * Q[5]) + x6 * (Q[6] + x * Q[7] + x2 * Q[8]);

		return remainder + multiplier * x * p / q;
	}

	/// @brief Generates a list of uniformly spaced values from `min` to `max`, inclusive, separated by `delta`.
	/// For example, `linear_space(1, 5, 1.0)` generates the vector `[1.0, 2.0, 3.0, 4.0, 5.0]`.
	/// @param[in] min Start of the list
	/// @param[in] max End of the list
	/// @param[in] delta Step size
	inline std::vector<double> linear_space(const double min, const double max, const double delta) {
		const std::size_t count = static_cast<std::size_t>(std::floor((max - min) / delta) + 1.0);

		std::vector<double> list(count);

		for (std::size_t i = 0; i < count; i++) {
			const double value = min + static_cast<double>(i) * delta;
			list[i] = value;
		}

		return list;
	}

	/// @brief Generates a list of logarithmically-spaced values from `min` to `max`, inclusive, separated by `delta` (in log space).
	/// For example, `log_space(1e1, 1e5, 1.0)` generates the vector `[1e1, 1e2, 1e3, 1e4, 1e5]`.
	/// @param[in] min Start of the list
	/// @param[in] max End of the list
	/// @param[in] delta Step in log-space
	inline std::vector<double> log_space(const double min, const double max, const double delta) {
		const double min_log = std::log10(min);
		const double max_log = std::log10(max);
		const std::size_t count = static_cast<std::size_t>(std::floor((max_log - min_log) / delta) + 1.0);

		std::vector<double> list(count);

		for (std::size_t i = 0; i < count; i++) {
			const double log10 = min_log + static_cast<double>(i) * delta;
			const double value = std::pow(10, log10);
			list[i] = value;
		}

		return list;
	}

	/// @brief Generates a list of `count` Chebyshev points from `min` to `max`, inclusive, in increasing order: the extrema of the
	/// Chebyshev polynomial of degree `count - 1`, mapped from [-1, 1] to [min, max]
	/// @param[in] min Start of the list
	/// @param[in] max End of the list
	/// @param[in] count Number of points (must be at least 2)
	inline std::vector<double> chebyshev_space(const double min, const double max, const std::size_t count) {
		std::vector<double> list(count);

		for (std::size_t i = 0; i < count; i++) {
			const double angle = std::numbers::pi * static_cast<double>(i) / static_cast<double>(count - 1);
			const double value = 0.5 * (min + max) - 0.5 * (max - min) * std::cos(angle);
			list[i] = value;
		}

		return list;
	}

	/// Checks whether the given value \p v is in the range [`l`, `h`] (inclusive)
	inline bool in_range(const auto v, const auto l, const auto h) noexcept {
		return !(v < l) && !(h < v);
	}
	/// @brief Converts the given boolean into a number (`true` ==> `1.0`, `false` ==> `0.0`)
	/// @details Since booleans are already represented by either a 0 or a 1, this function is just a cast
	inline double boole(const bool v) noexcept {
		return static_cast<double>(v);
	}

}

}

#endif