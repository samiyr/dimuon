#ifndef DM_INTERPOLATION_LINEAR_INTERPOLATION_BASIS_H
#define DM_INTERPOLATION_LINEAR_INTERPOLATION_BASIS_H

#include <algorithm>
#include <cstddef>
#include <functional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "dimuon/interpolation/interpolation_basis.hpp"
#include "dimuon/utility/shims.hpp" // IWYU pragma: keep
#include "dimuon/utility/collections.hpp"

namespace dm {

/// The kernels of piecewise-linear interpolation on a list of knots, with one triangular ("hat") kernel for each knot.
/// A single kernel `k` is one at `knots[k]`, falls linearly to zero at `knots[k - 1]` and `knots[k + 1]`, and is zero outside them.
/// The kernels of the first and the last knot are only half of a full triangular kernel: they reach one at the edge, but only fall off on one side.
/// At any point at most two kernels are non-zero, and between the first and the last knot the kernels sum to one,
/// so that the sum of `f(knots[k])` times kernel `k` is the piecewise-linear interpolation of `f` over the whole range of the knots.
/// The kernels use the mapped variable `t = mapping(x)`. The default is the identity map.
/// @tparam mapping_t Type of the mapping applied via `t = mapping(x)`, defaults to the identity map
template <typename mapping_t = std::identity>
class linear_interpolation_basis {
public:
	/// The kernels are piecewise-linear in `t`
	static constexpr std::size_t degree = 1;
	/// The type of the mapping `t = mapping(x)`
	using mapping_type = mapping_t;
	/// At most two kernels can be non-zero at any given point (see illustration of `linear_interpolation_basis::nonzero_at`)
	static constexpr size_t maximum_kernels_per_point = 2;
	/// Triangular kernel, which is one at `middle_knot` and falls linearly (in `t`) to zero at the knots on either side of it, and zero outside them.
	/// A half kernel has no knot on one side. One of the edge knots (`left_knot` or `right_knot`) is equal to `middle_knot`.
	struct kernel {
		kernel_index index;

		/// Constructs the kernel from its \p knots_ (left, middle and right) and their \p nodes_, the knots mapped to `t`
		explicit kernel(const kernel_index index_, const std::span<const double> knots_, const std::span<const double> nodes_) noexcept
		: index{index_},
		left_knot{knots_[0]},
		middle_knot{knots_[1]},
		right_knot{knots_[2]},
		left_node{nodes_[0]},
		middle_node{nodes_[1]},
		right_node{nodes_[2]} {}

		/// Returns the value of the kernel at \p x, which is mapped to \p t
		double operator()(const double x, const double t) const noexcept {
			if (!contains(x)) { return 0.0; }
			if (x < middle_knot) { return (t - left_node) / (middle_node - left_node); }
			if (x > middle_knot) { return 1.0 - (t - middle_node) / (right_node - middle_node); }
			return 1.0;
		}

		/// Returns whether the kernel is non-zero at \p x, which it is strictly between its outer knots and at its middle knot (an outer knot of a half kernel)
		bool contains(const double x) const noexcept { return (left_knot < x && x < right_knot) || x == middle_knot; }

	private:
		double left_knot;
		double middle_knot;
		double right_knot;

		double left_node;
		double middle_node;
		double right_node;
	};

	/// @brief Constructs the interpolation basis from the list of knots
	/// @param[in] knots_ The knots, in increasing order. There must be at least two.
	/// @param[in] mapping_ The mapping of the kernel variable (`t = mapping(x)`), which must be finite and strictly monotonic at the knots
	/// @param[in] index_offset Starting index for the numbering of kernel indices. Used by `composite_interpolation_basis`. Defaults to 0.
	explicit linear_interpolation_basis(std::vector<double> knots_, mapping_t mapping_ = {}, const kernel_index index_offset = 0) :
		knots{std::move(knots_)},
		mapping{std::move(mapping_)}
	{
		const std::vector<double> nodes = knots | std::views::transform(mapping) | std::ranges::to<std::vector>();
		check_interpolation_knots(knots, nodes);

		kernels = construct_kernels(knots, nodes, index_offset);
	}

	/// Check whether \p x is in the knot range
	bool in_range(const double x) const noexcept {
		return knots.front() <= x && x <= knots.back();
	}

	/// Returns the number of kernels
	std::size_t size() const noexcept { return kernels.size(); }

	/// Returns the kernels that are non-zero at \p x (at most two), as a range
	///																									///
	///                               │                                                   				///
	///  \     /\     /\     /\     /\│    /\     /\     /\     /\     /\     /\     /    				///
	///   \   /  \   /  \   /  \   /  │   /  \   /  \   /  \   /  \   /  \   /  \   /     				///
	///    \ /    \ /    \ /    \ /   │\ /    \ /    \ /    \ /    \ /    \ /    \ /      				///
	///     /      \      /      \    │ /      \      /      \      /      \      /       				///
	///    / \    / \    / \    / \   │/ \    / \    / \    / \    / \    / \    / \      				///
	///   /   \  /   \  /   \  /   \  │   \  /   \  /   \  /   \  /   \  /   \  /   \     				///
	///  /     \/     \/     \/     \/│    \/     \/     \/     \/     \/     \/     \    				///
	/// ──────────────────────────────┼──────────────────────────────────────────────────►				///
	///                               │                                                   				///
	///																									///
	auto nonzero_at(const double x) const noexcept {
		// Only the two kernels that peak at the ends of the knot interval that contains x can be non-zero
		using diff_t = typename decltype(knots)::difference_type;

		const diff_t interval = collections::iterator_to_index(knots, std::ranges::upper_bound(knots, x));
		const diff_t start_idx = std::max(interval - 1, static_cast<diff_t>(0));
		const diff_t end_idx = interval + 1;
		
		return kernels | shims::views::slice(start_idx, end_idx) | std::views::filter([x](const kernel &k) { return k.contains(x); });
	}

	/// Returns the values of the kernels that are non-zero at \p x (via `nonzero_at`) as a range of (index, value) pairs
	auto values_at(const double x) const noexcept {
		const double t = mapping(x);
		return nonzero_at(x) | std::views::transform([x, t](const kernel &k) { return std::pair{k.index, k(x, t)}; });
	}

private:
	/// Returns the kernels of \p knots_ and their \p nodes_, numbered from \p index_offset: kernel `k` is made of `knots_[k]` and the knots on either side of it
	static std::vector<kernel> construct_kernels(const std::vector<double> &knots_, const std::vector<double> &nodes_, const kernel_index index_offset) {
		// The first and the last knot have no knot on one side. They are repeated there, which makes their kernels halves (see `kernel`).
		const auto pad = [](const std::vector<double> &points) {
			std::vector<double> padded_points{points.front()};
			padded_points.append_range(points);
			padded_points.push_back(points.back());
			return padded_points;
		};
		const std::vector<double> padded_knots = pad(knots_);
		const std::vector<double> padded_nodes = pad(nodes_);

		return std::views::zip_transform(
			[](const kernel_index index, const std::span<const double> knot_window, const std::span<const double> node_window) {
				return kernel{index, knot_window, node_window};
			},
			std::views::iota(index_offset),
			padded_knots | std::views::slide(3),
			padded_nodes | std::views::slide(3)
		) | std::ranges::to<std::vector>();
	}

	/// The knots, in increasing order
	std::vector<double> knots;
	/// The mapping of the kernel variable
	mapping_t mapping;
	/// The kernels, one for each knot
	std::vector<kernel> kernels;
};

}

#endif
