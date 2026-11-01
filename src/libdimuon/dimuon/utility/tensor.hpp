#ifndef DM_UTILITY_TENSOR_H
#define DM_UTILITY_TENSOR_H

#include <cassert>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "dimuon/utility/shims.hpp"

namespace dm {

/// @brief A tensor of rank \p rank with elements of type \p value_t, indexed as `tensor[i, j, ...]` with one index for each dimension.
/// The elements are stored in a container of type \p storage_t in row-major order (`shims::layout_right`), and the tensor is a view of them (`shims::mdspan`, see `view`).
/// @details This is a simple version of `std::mdarray`
/// @tparam storage_t Type of the container of the elements: `std::vector` (the default), or e.g. `shims::inplace_vector` for a tensor that allocates nothing
template <typename value_t, std::size_t rank, typename storage_t = std::vector<value_t>>
class tensor {
public:
	/// The size of each dimension
	using extents_t = shims::dextents<std::size_t, rank>;

	tensor() = default;

	/// Constructs a tensor with the given extents and every element set to \p value
	explicit tensor(const extents_t &extents, const value_t &value) : layout{extents}, elements(layout.required_span_size(), value) {}

	/// Constructs a tensor with the given extents and elements in row-major order.
	/// The dimensions must match (i.e. `elements.size()` is the product of the extents), otherwise future usage of this class may result in undefined behaviour.
	explicit tensor(const extents_t &extents, storage_t elements_) noexcept(std::is_nothrow_move_constructible_v<storage_t>) : layout{extents}, elements{std::move(elements_)} {
		assert(elements.size() == layout.required_span_size());
	}

	/// Returns the size of \p dimension
	std::size_t extent(const std::size_t dimension) const noexcept { return layout.extents().extent(dimension); }

	/// Returns a reference to the element with the given indices, one for each dimension
	auto &operator[](this auto &self, const std::convertible_to<std::size_t> auto... indices) noexcept requires (sizeof...(indices) == rank) {
		return self.view()[indices...];
	}

	/// The elements in row-major order
	auto begin(this auto &self) noexcept { return self.elements.begin(); }
	/// The elements in row-major order
	auto end(this auto &self) noexcept { return self.elements.end(); }

	/// Returns the tensor as a view of its elements. The view is made here and is not a member, since in a copy of the tensor,
	/// a member view would still be a view of the elements of the original.
	auto view(this auto &self) noexcept { return shims::mdspan{self.elements.data(), self.layout}; }

private:
	/// The row-major layout of the elements, which maps the indices of an element to its position in `elements`
	shims::layout_right::mapping<extents_t> layout{};
	/// The elements in row-major order
	storage_t elements;
};

}

#endif
