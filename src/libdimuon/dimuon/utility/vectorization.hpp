#ifndef DM_UTILITY_VECTORIZATION_H
#define DM_UTILITY_VECTORIZATION_H

#include <concepts>
#include <cstddef>
#include <experimental/simd> // IWYU pragma: keep (clangd does not map std::experimental::simd to it)
#include <tuple>
#include <array>
#include <algorithm>

namespace dm {

/// Supporting functionality for SIMD vectorization
namespace vectorization {
	/// Vector of \p lanes doubles, with element-wise arithmetic (see `std::experimental::simd`)
	template <std::size_t lanes>
	using vec = std::experimental::fixed_size_simd<double, lanes>;

	/// A concept of a vector of doubles such that
	/// - the number of lanes is known at compile time
	/// - the value in lane `lane` is given by `values[lane]`.
	/// Satisfied by a SIMD vector (e.g. `vec`), or an array such as `std::array` of doubles.
	template <typename vector_t>
	concept fixed_size_vector = (std::experimental::is_simd_v<vector_t> || requires { std::tuple_size<vector_t>::value; })
		&& requires(const vector_t &values, const std::size_t lane) { { values[lane] } -> std::convertible_to<double>; };

	static_assert(fixed_size_vector<vec<1>>, "vectorization::vec does not satisfy the concept fixed_size_vector");
	static_assert(fixed_size_vector<std::array<double, 1>>, "std::array does not satisfy the concept fixed_size_vector");

	/// Gives the number of lanes of the fixed-size vector type \p vector_t
	template <fixed_size_vector vector_t>
	inline constexpr std::size_t vec_size_v = [] {
		if constexpr (std::experimental::is_simd_v<vector_t>) { 
			return vector_t::size();
		} else {
			return std::tuple_size_v<vector_t>;
		}
	}();

	/// @brief Computes the Kronecker product `a * b` of \p a and \p b, i.e. the products in the order of the pairs `(i, j)`, with `j` varying the fastest
	/// @return A vector of all products `a[i] * b[j]`, where `a[i] * b[j]` is in lane `i * vec_size_v<vec_b_t> + j`
	template <fixed_size_vector vec_a_t, fixed_size_vector vec_b_t>
	vec<vec_size_v<vec_a_t> * vec_size_v<vec_b_t>> kronecker_product(const vec_a_t &a, const vec_b_t &b) noexcept {
		constexpr std::size_t size_b = vec_size_v<vec_b_t>;
		return vec<vec_size_v<vec_a_t> * size_b>([&](const std::size_t lane) { return a[lane / size_b] * b[lane % size_b]; });
	}

	/// Constructs a \p lanes -dimensional standard basis vector, 
	/// i.e. the value of lane \p index is one and the others are zero
	template <std::size_t lanes>
	vec<lanes> unit_vector(const std::size_t index) noexcept {
		return vec<lanes>([index](const std::size_t lane) { return lane == index ? 1.0 : 0.0; });
	}

	/// Constructs the \p lanes -dimensional identity vector with respect to the Hadamard elementwise product,
	/// i.e. the value of every lane is one
	template <std::size_t lanes>
	vec<lanes> ones() noexcept {
		return vec<lanes>(1.0);
	}

	/// @brief Sums the lanes of the vector
	/// @details If `vector_t` is `vec` (i.e. a SIMD vector), then the summation order is unspecified
	template <fixed_size_vector vector_t>
	double sum_lanes(const vector_t &values) noexcept {
		if constexpr (std::experimental::is_simd_v<vector_t>) { 
			return std::experimental::reduce(values);
		} else {
			return std::ranges::fold_left(values, 0.0, std::plus{});
		}
	}

	/// Computes the dot product of two vectors
	template <fixed_size_vector vec_a_t, fixed_size_vector vec_b_t>
		requires (vec_size_v<vec_a_t> == vec_size_v<vec_b_t>)
	double dot_product(const vec_a_t &a, const vec_b_t &b) noexcept {
		return sum_lanes(vec<vec_size_v<vec_a_t>>([&](const std::size_t lane) { return a[lane] * b[lane]; }));
	}

	/// @brief Calculates the contraction defined in eq. (14) of `numerical_evaluation.md`
	/// @return A vector with the sum of `a[m] * b[d * m + i]` over `m` in lane `i`
	template <fixed_size_vector vec_a_t, fixed_size_vector vec_b_t>
		requires (vec_size_v<vec_b_t> % vec_size_v<vec_a_t> == 0)
	vec<vec_size_v<vec_b_t> / vec_size_v<vec_a_t>> contract(const vec_a_t &a, const vec_b_t &b) noexcept {
		constexpr std::size_t d_1 = vec_size_v<vec_a_t>;
		constexpr std::size_t d_2 = vec_size_v<vec_b_t>;
		constexpr std::size_t d = d_2 / d_1;

		return vec<d>([&](const std::size_t i) {
			double sum = 0.0;
			for (std::size_t m = 0; m < d_1; m++) { sum += a[m] * b[d * m + i]; }
			return sum;
		});
	}
}

}

#endif
