#ifndef DM_INTERPOLATION_BLOCK_INTERPOLATION_BASIS_H
#define DM_INTERPOLATION_BLOCK_INTERPOLATION_BASIS_H

#include <algorithm>
#include <cstddef>
#include <format>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dimuon/interpolation/interpolation_basis.hpp"
#include "dimuon/utility/collections.hpp"

namespace dm {

/// A block-based interpolation basis. The entire interpolation domain is grouped into overlapping blocks,
/// and interpolation takes place within blocks. Each block consists of `degree + 1` knots.
/// The actual interpolation is handled by the supplied block basis.
/// @tparam block_basis_t Type of interpolation basis used within each block
/// @details The approach here follows EKO (https://eko.readthedocs.io/en/latest/theory/Interpolation.html)
template <interpolation_basis block_basis_t>
class block_interpolation_basis {
public:
	/// The degree of the blocks
	static constexpr std::size_t degree = block_basis_t::degree;
	/// The kernel variable mapping of the blocks
	using mapping_type = typename block_basis_t::mapping_type;

	static_assert(degree != std::dynamic_extent, "Block interpolation needs a fixed-degree basis");
	static_assert(degree >= 1, "Block interpolation needs a basis of degree one or higher");

	/// Only the kernels of a single block can be non-zero at any given point
	static constexpr std::size_t maximum_kernels_per_point = block_basis_t::maximum_kernels_per_point;

	/// @brief Constructs a block-based interpolation from the list of knots
	/// @param[in] knots_ The knots, in increasing order. There must be at least `degree + 1`.
	/// @param[in] mapping_ The mapping of the kernel variable (`t = mapping(x)`), which must be finite and strictly monotonic at the knots
	explicit block_interpolation_basis(std::vector<double> knots_, mapping_type mapping_ = {}) :
		knots{std::move(knots_)},
		mapping{std::move(mapping_)},
		blocks{construct_blocks(knots, mapping)}
	{
		const std::vector<double> nodes = knots | std::views::transform(mapping) | std::ranges::to<std::vector>();
		check_interpolation_knots(knots, nodes);
		if (knots.size() < degree + 1) {
			throw std::invalid_argument(std::format("Block interpolation of degree {} needs at least {} knots, got {}", degree, degree + 1, knots.size()));
		}
	}

	/// Returns the number of kernels
	std::size_t size() const noexcept { return knots.size(); }

	/// Check whether \p x is in the knot range
	bool in_range(const double x) const noexcept {
		return knots.front() <= x && x <= knots.back();
	}

	/// Returns the kernels that are non-zero at \p x, as a range.
	/// The kernels are those of the active block, as determined by \p x.
	auto nonzero_at(const double x) const {
		return block_at(x).nonzero_at(x);
	}

	/// Returns the values of the kernels that are non-zero at \p x (via `nonzero_at`) as a range of (index, value) pairs
	auto values_at(const double x) const {
		return block_at(x).values_at(x);
	}

private:
	/// Returns the active block at \p x.
	/// Any given \p x can have multiple active block candidates, so the block is chosen (as in EKO) as the middle block.
	/// When there are an even number of blocks to choose from (i.e. no genuine middle block), the "upper" one is chosen.
	/// Near the edges, the choice can be clamped to one of the edges instead.
	const block_basis_t &block_at(const double x) const noexcept {
		if (!in_range(x)) { return blocks.front(); }

		using diff_t = typename decltype(knots)::difference_type;

		// The index of the knot before the first knot that contains `x`
		const diff_t knot_index = std::clamp(collections::iterator_to_index(knots, std::ranges::lower_bound(knots, x)) - 1, diff_t{0}, std::ssize(knots) - 2);
		// The index of the block where `x` is most central (when possible)
		const diff_t block_index = std::clamp(knot_index - static_cast<diff_t>((degree - 1) / 2), diff_t{0}, std::ssize(blocks) - 1);

		return blocks[static_cast<std::size_t>(block_index)];
	}

	/// Returns the blocks of \p knots_, with the correct index offsets applied
	static std::vector<block_basis_t> construct_blocks(const std::vector<double> &knots_, const mapping_type &mapping_) {
		return std::views::zip_transform(
			[&mapping_](const kernel_index first_kernel_index, const std::span<const double> window) -> block_basis_t {
				return block_basis_t{window | std::ranges::to<std::vector>(), mapping_, first_kernel_index};
			},
			std::views::iota(kernel_index{0}),
			knots_ | std::views::slide(degree + 1)
		) | std::ranges::to<std::vector>();
	}

	/// The knots, in increasing order
	std::vector<double> knots;
	/// The mapping of the kernel variable
	mapping_type mapping;
	/// The blocks
	std::vector<block_basis_t> blocks;
};

}

#endif
