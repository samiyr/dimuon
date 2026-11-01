#ifndef DM_UTILITY_COLLECTIONS_H
#define DM_UTILITY_COLLECTIONS_H

#include <concepts>
#include <vector>
#include <ranges>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <cstddef>
#include <type_traits>

namespace dm {

namespace collections {
	template <typename value_t>
	inline std::vector<value_t> join(const std::vector<value_t> &a, const std::vector<value_t> &b) {
		std::vector<value_t> output;
		output.reserve(a.size() + b.size());
		output.insert(output.end(), a.cbegin(), a.cend());
		output.insert(output.end(), b.cbegin(), b.cend());
		return output;
	}

	template <typename value_t, typename... args_t>
	inline std::vector<value_t> join(const std::vector<value_t> &a, const args_t&... b) {
		return join(a, join(b...));
	}

	template <typename value_t> 
	inline void hash_combine(std::size_t& seed, const value_t &v) noexcept {
		std::hash<value_t> hasher;
		const std::size_t kMul = 0x9ddfea08eb382d69ULL;
		std::size_t a = (hasher(v) ^ seed) * kMul;
		a ^= (a >> 47);
		std::size_t b = (seed ^ a) * kMul;
		b ^= (b >> 47);
		seed = b * kMul;
	}

	/// Checks whether the **sorted** range \p r contains only unique values
	inline bool is_unique(const std::ranges::range auto &r) noexcept {
		return std::ranges::adjacent_find(r) == std::ranges::end(r);
	}

	template <typename R, typename V>
	concept range_of_value = std::ranges::range<R> && std::same_as<std::ranges::range_value_t<R>, V>;

	/// A range whose elements are convertible to `V`
	template <typename R, typename V>
	concept range_convertible_to = std::ranges::range<R> && std::convertible_to<std::ranges::range_reference_t<R>, V>;

	inline auto iterator_to_index(const auto &cont, const auto &it) noexcept {
		return std::distance(std::cbegin(cont), it);
	}

	// Element-wise arithmetic of `std::array`. `std::array` is not in this namespace, so the operators are not found by argument-dependent lookup:
	// an expression finds them only where they have been declared, e.g. with `using collections::operator+;`.

	/// Adds \p rhs to \p lhs element by element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator+=(std::array<value_t, size> &lhs, const std::array<value_t, size> &rhs) noexcept {
		for (auto &&[left, right] : std::views::zip(lhs, rhs)) { left += right; }
		return lhs;
	}

	/// Subtracts \p rhs from \p lhs element by element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator-=(std::array<value_t, size> &lhs, const std::array<value_t, size> &rhs) noexcept {
		for (auto &&[left, right] : std::views::zip(lhs, rhs)) { left -= right; }
		return lhs;
	}

	/// Multiplies \p lhs by \p rhs element by element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator*=(std::array<value_t, size> &lhs, const std::array<value_t, size> &rhs) noexcept {
		for (auto &&[left, right] : std::views::zip(lhs, rhs)) { left *= right; }
		return lhs;
	}

	/// Divides \p lhs by \p rhs element by element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator/=(std::array<value_t, size> &lhs, const std::array<value_t, size> &rhs) noexcept {
		for (auto &&[left, right] : std::views::zip(lhs, rhs)) { left /= right; }
		return lhs;
	}

	/// Adds \p rhs, or its conversion to `value_t`, to each element of \p lhs
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator+=(std::array<value_t, size> &lhs, const std::type_identity_t<value_t> &rhs) noexcept {
		for (value_t &value : lhs) { value += rhs; }
		return lhs;
	}

	/// Subtracts \p rhs, or its conversion to `value_t`, from each element of \p lhs
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator-=(std::array<value_t, size> &lhs, const std::type_identity_t<value_t> &rhs) noexcept {
		for (value_t &value : lhs) { value -= rhs; }
		return lhs;
	}

	/// Multiplies each element of \p lhs by \p rhs, or by its conversion to `value_t`
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator*=(std::array<value_t, size> &lhs, const std::type_identity_t<value_t> &rhs) noexcept {
		for (value_t &value : lhs) { value *= rhs; }
		return lhs;
	}

	/// Divides each element of \p lhs by \p rhs, or by its conversion to `value_t`
	template <typename value_t, std::size_t size>
	std::array<value_t, size> &operator/=(std::array<value_t, size> &lhs, const std::type_identity_t<value_t> &rhs) noexcept {
		for (value_t &value : lhs) { value /= rhs; }
		return lhs;
	}

	/// Element-wise sum
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator+(std::array<value_t, size> lhs, const std::array<value_t, size> &rhs) noexcept { lhs += rhs; return lhs; }
	/// Element-wise difference
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator-(std::array<value_t, size> lhs, const std::array<value_t, size> &rhs) noexcept { lhs -= rhs; return lhs; }
	/// Element-wise product
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator*(std::array<value_t, size> lhs, const std::array<value_t, size> &rhs) noexcept { lhs *= rhs; return lhs; }
	/// Element-wise quotient
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator/(std::array<value_t, size> lhs, const std::array<value_t, size> &rhs) noexcept { lhs /= rhs; return lhs; }

	/// Sum of each element and \p rhs
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator+(std::array<value_t, size> lhs, const std::type_identity_t<value_t> &rhs) noexcept { lhs += rhs; return lhs; }
	/// Sum of \p lhs and each element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator+(const std::type_identity_t<value_t> &lhs, std::array<value_t, size> rhs) noexcept { rhs += lhs; return rhs; }
	/// Product of each element and \p rhs
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator*(std::array<value_t, size> lhs, const std::type_identity_t<value_t> &rhs) noexcept { lhs *= rhs; return lhs; }
	/// Product of \p lhs and each element
	template <typename value_t, std::size_t size>
	std::array<value_t, size> operator*(const std::type_identity_t<value_t> &lhs, std::array<value_t, size> rhs) noexcept { rhs *= lhs; return rhs; }
}

}

#endif