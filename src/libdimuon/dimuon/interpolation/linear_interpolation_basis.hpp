#ifndef DM_INTERPOLATION_LINEAR_INTERPOLATION_BASIS_H
#define DM_INTERPOLATION_LINEAR_INTERPOLATION_BASIS_H

#include <algorithm>
#include <cstddef>
#include <format>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "dimuon/utility/collections.hpp"
#include "dimuon/utility/shims.hpp" // IWYU pragma: keep (provides shims::views::slice, which clangd attributes to slice_view.hpp)

namespace dm {

/// Index of an interpolation kernel (see `linear_interpolation_basis`)
using kernel_index = std::size_t;

/// The kernels of piecewise-linear interpolation on a list of knots, with one triangular ("hat") kernel for each knot.
/// A single kernel `k` is one at `knots[k]`, falls linearly to zero at `knots[k - 1]` and `knots[k + 1]`, and is zero outside them.
/// The kernels of the first and the last knot are only half of a full triangular kernel: they reach one at the edge, but only fall off on one side.
/// At any point at most two kernels are non-zero, and between the first and the last knot the kernels sum to one,
/// so that the sum of `f(knots[k])` times kernel `k` is the piecewise-linear interpolation of `f` over the whole range of the knots.
class linear_interpolation_basis {
public:
	static constexpr size_t maximum_kernels_per_point = 2;
	/// Triangular kernel, which is one at `middle_knot` and falls linearly to zero at the knots on either side of it, and zero outside them.
	/// A half kernel has no knot on one side. One of the edge knots (`left_knot` or `right_knot`) is equal to `middle_knot`.
	struct kernel {
		kernel_index index;

		explicit kernel(const kernel_index index_, const double left_knot_, const double middle_knot_, const double right_knot_) noexcept
		: index{index_},
		left_knot{left_knot_},
		middle_knot{middle_knot_},
		right_knot{right_knot_} {}

		/// Returns the value of the kernel at \p x
		double operator()(const double x) const noexcept {
			if (!contains(x)) { return 0.0; }
			if (x < middle_knot) { return (x - left_knot) / (middle_knot - left_knot); }
			if (x > middle_knot) { return 1.0 - (x - middle_knot) / (right_knot - middle_knot); }
			return 1.0;
		}

		/// Returns whether the kernel is non-zero at \p x, which it is strictly between its outer knots and at its middle knot (an outer knot of a half kernel)
		bool contains(const double x) const noexcept { return (left_knot < x && x < right_knot) || x == middle_knot; }

	private:
		double left_knot;
		double middle_knot;
		double right_knot;
	};

	/// @brief Constructs the interpolation basis from the list of knots
	/// @param[in] knots_ The knots, in increasing order. There must be at least two.
	explicit linear_interpolation_basis(std::vector<double> knots_) : knots{std::move(knots_)} {
		if (knots.size() < 2) {
			throw std::invalid_argument(std::format("At least two interpolation knots must be supplied, got {}", knots.size()));
		}
		if (!(std::ranges::is_sorted(knots) && knots.front() < knots.back())) {
			throw std::invalid_argument("Interpolation knots must be in increasing order");
		}
		if (!collections::is_unique(knots)) {
			throw std::invalid_argument("Interpolation knots must be unique");
		}

		kernels = construct_kernels(knots);
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
		using diff_t = decltype(knots)::difference_type;

		const diff_t interval = std::ranges::upper_bound(knots, x) - knots.begin();
		const diff_t start_idx = std::max(interval - 1, static_cast<diff_t>(0));
		const diff_t end_idx = interval + 1;
		
		return kernels | shims::views::slice(start_idx, end_idx) | std::views::filter([x](const kernel &k) { return k.contains(x); });
	}

private:
	/// Returns the kernels of \p knots_: kernel `k` is made of `knots_[k]` and the knots on either side of it
	static std::vector<kernel> construct_kernels(const std::vector<double> &knots_) {
		// The first and the last knot have no knot on one side. They are repeated there, which makes their kernels halves (see `kernel`).
		std::vector<double> padded_knots{knots_.front()};
		padded_knots.append_range(knots_);
		padded_knots.push_back(knots_.back());

		return std::views::zip_transform(
			[](const kernel_index index, const std::span<const double> window) { return kernel{index, window[0], window[1], window[2]}; },
			std::views::iota(kernel_index{0}),
			padded_knots | std::views::slide(3)
		) | std::ranges::to<std::vector>();
	}

	std::vector<double> knots;
	std::vector<kernel> kernels;
};

}

#endif
