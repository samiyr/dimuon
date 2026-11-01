#ifndef DM_UTILITY_SPARSE_VECTOR_H
#define DM_UTILITY_SPARSE_VECTOR_H

#include <cstddef>
#include <flat_map>
#include <functional>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "dimuon/utility/shims.hpp"

namespace dm {

/// @brief A simple sparse vector, which stores only the non-zero elements (non-zero being a value `!= value_t{}`).
/// Iteration gives the stored elements as (index, value) pairs, in increasing order of the index.
/// @details The sparse vector is implemented with `std::flat_map`, which inserts `value_t{}` whenever a non-existing key 
/// is accessed. This means that accessing a non-existing index modifies the sparse vector by inserting a `value_t{}` at that
/// index. It may therefore also contain zero elements (i.e. `value_t{}`).
/// @tparam value_t Type of the elements
/// @tparam capacity The largest number of elements that can be stored, or `std::dynamic_extent` (the default) for no limit.
/// With a limit, the elements are stored in place (see `shims::inplace_vector`), so the sparse vector allocates nothing.
template <typename value_t, std::size_t capacity = std::dynamic_extent>
class sparse_vector {
public:
	/// Sets all elements to `value_t{}`
	void clear() noexcept { elements.clear(); }

	/// @brief Returns whether no elements are stored
	/// @details Due to the behaviour of `std::flat_map`, an empty sparse vector will return `value_t{}` at any index
	bool empty() const noexcept { return elements.empty(); }

	/// Returns a reference to the element at \p index if it has previously been assigned, 
	/// otherwise `value_t{}` is assigned to the index and a reference to it is returned
	value_t &operator[](const std::size_t index) { return elements[index]; }

	/// The stored elements as (index, value) pairs, in increasing order of the index
	auto begin() const noexcept { return elements.begin(); }
	/// The stored elements as (index, value) pairs, in increasing order of the index
	auto end() const noexcept { return elements.end(); }

	/// Adds \p rhs element by element. The sum stores the elements that either of the two stores.
	sparse_vector &operator+=(const sparse_vector &rhs) {
		for (const auto &[index, value] : rhs) { elements[index] += value; }
		return *this;
	}

	/// Subtracts \p rhs element by element. The difference stores the elements that either of the two stores.
	sparse_vector &operator-=(const sparse_vector &rhs) {
		for (const auto &[index, value] : rhs) { elements[index] -= value; }
		return *this;
	}

	/// Multiplies by \p rhs element by element. Each stored element is multiplied by the element of \p rhs at the same index,
	/// which is `value_t{}` if \p rhs stores none, so the product stores the same elements as before.
	sparse_vector &operator*=(const sparse_vector &rhs) noexcept {
		for (auto &&[index, value] : elements) { value *= rhs.element_or_zero(index); }
		return *this;
	}

	/// Divides by \p rhs element by element. Each stored element is divided by the element of \p rhs at the same index,
	/// which is `value_t{}` if \p rhs stores none, so the quotient stores the same elements as before.
	/// The elements that are not stored stay zero, also where the element of \p rhs is zero.
	sparse_vector &operator/=(const sparse_vector &rhs) noexcept {
		for (auto &&[index, value] : elements) { value /= rhs.element_or_zero(index); }
		return *this;
	}

	/// Adds \p rhs to each stored element. The elements that are not stored stay zero,
	/// so the sum depends on which elements are stored and not only on their values.
	sparse_vector &operator+=(const value_t &rhs) noexcept {
		for (value_t &value : elements | std::views::values) { value += rhs; }
		return *this;
	}

	/// Multiplies each element by \p rhs
	sparse_vector &operator*=(const value_t &rhs) noexcept {
		for (value_t &value : elements | std::views::values) { value *= rhs; }
		return *this;
	}

	/// Element-wise sum (see `operator+=`)
	friend sparse_vector operator+(sparse_vector lhs, const sparse_vector &rhs) { lhs += rhs; return lhs; }
	/// Element-wise difference (see `operator-=`)
	friend sparse_vector operator-(sparse_vector lhs, const sparse_vector &rhs) { lhs -= rhs; return lhs; }
	/// Element-wise product (see `operator*=`)
	friend sparse_vector operator*(sparse_vector lhs, const sparse_vector &rhs) noexcept { lhs *= rhs; return lhs; }
	/// Element-wise quotient (see `operator/=`)
	friend sparse_vector operator/(sparse_vector lhs, const sparse_vector &rhs) noexcept { lhs /= rhs; return lhs; }

	/// Sum of each stored element and \p rhs (see `operator+=`)
	friend sparse_vector operator+(sparse_vector lhs, const value_t &rhs) noexcept { lhs += rhs; return lhs; }
	/// Sum of \p lhs and each stored element (see `operator+=`)
	friend sparse_vector operator+(const value_t &lhs, sparse_vector rhs) noexcept { rhs += lhs; return rhs; }
	/// Product of each element and \p rhs
	friend sparse_vector operator*(sparse_vector lhs, const value_t &rhs) noexcept { lhs *= rhs; return lhs; }
	/// Product of \p lhs and each element
	friend sparse_vector operator*(const value_t &lhs, sparse_vector rhs) noexcept { rhs *= lhs; return rhs; }

private:
	/// Container of the indices or the values of the stored elements: `std::vector`, or with a limited capacity, an in-place vector
	template <typename element_t>
	using container_t = std::conditional_t<
		capacity == std::dynamic_extent, std::vector<element_t>, shims::inplace_vector<element_t, capacity>
	>;
	
	using key_container_t = container_t<std::size_t>;
	using mapped_container_t = container_t<value_t>;

	/// Underlying map, which stores the values by index
	std::flat_map<std::size_t, value_t, std::less<std::size_t>, key_container_t, mapped_container_t> elements;

	/// Returns the element at \p index, or `value_t{}` if none is stored. Unlike `operator[]`, this stores nothing.
	value_t element_or_zero(const std::size_t index) const noexcept {
		const auto element = elements.find(index);
		return element != elements.end() ? element->second : value_t{};
	}
};

/// Returns the vector whose only non-zero element is \p value at \p index, or the zero vector if \p value is zero,
/// as a range of (index, value) pairs like a `sparse_vector`. Unlike a `sparse_vector`, it holds the element itself, so it allocates nothing.
template <typename value_t>
auto single_element_sparse_vector(const std::size_t index, const value_t &value) noexcept {
	return std::views::repeat(std::pair{index, value}, value == value_t{} ? 0 : 1);
}

}

#endif
