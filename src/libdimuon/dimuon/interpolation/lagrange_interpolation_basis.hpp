#ifndef DM_INTERPOLATION_LAGRANGE_INTERPOLATION_BASIS_H
#define DM_INTERPOLATION_LAGRANGE_INTERPOLATION_BASIS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <functional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dimuon/interpolation/linear_interpolation_basis.hpp"
#include "dimuon/utility/collections.hpp"
#include "dimuon/utility/shims.hpp" // IWYU pragma: keep

namespace dm {

/// The kernels of global Lagrange interpolation on a list of knots, with one Lagrange kernel for each knot.
/// Kernel `k` is the Lagrange polynomial of `knots[k]` in the mapped momentum fraction `t = mapping(x)`.
/// This is a the polynomial of degree `size() - 1` in `t` which is one at `knots[k]` and zero at all the other knots.
/// Between knots, the value of a Lagrange polynomial is non-zero.
/// As with `linear_interpolation_basis`, the kernels are zero outside the first and the last knot.
/// Within the range of knots the kernels sum to one, so that the sum of `f(knots[k])` times kernel `k` is the 
/// polynomial interpolation of `f` (in `t`) over all the knots. 
/// Unlike the piecewise-linear kernels, every kernel is non-zero almost everywhere within the knot range.
/// Polynomial interpolation through many points is accurate and stable only for points that cluster towards the ends of the range, such as Chebyshev points in `t`.
/// Otherwise one runs into Runge's phenomenon.
/// @tparam mapping_t Type of the mapping of the momentum fraction to the variable `t` in which the kernels are polynomials (defaults to the identity map)
template <typename mapping_t = std::identity>
class lagrange_interpolation_basis {
public:
	/// All the kernels can be non-zero at a given point
	static constexpr std::size_t maximum_kernels_per_point = std::dynamic_extent;

	/// Lagrange kernel, the polynomial in the mapped momentum fraction that is one at its node and zero at all the other nodes.
	/// Its value is given by the numerically stable barycentric formula, `(w_k / (t - t_k)) / sum_j (w_j / (t - t_j))`, where `w_j` are the barycentric weights.
	struct kernel {
		kernel_index index;

		explicit kernel(const kernel_index index_, const double node_, const double weight_) noexcept : index{index_}, node{node_}, weight{weight_} {}

		/// Returns the value of the kernel at the mapped momentum fraction \p t given the sum of the barycentric weights.
		/// At the node (i.e. if `t == node`), return `1.0`.
		double operator()(const double t, const double weight_sum) const noexcept {
			if (t == node) { return 1.0; }
			return weight / (t - node) / weight_sum;
		}

	private:
		double node;
		double weight;
	};

	/// @brief Constructs the interpolation basis from the list of knots
	/// @param[in] knots_ The knots, in increasing order. There must be at least two.
	/// @param[in] mapping_ The mapping of the momentum fraction to the variable in which the kernels are polynomials, which must be finite and strictly
	/// monotonic at the knots
	explicit lagrange_interpolation_basis(std::vector<double> knots_, mapping_t mapping_ = {}) :
		knots{std::move(knots_)},
		mapping{std::move(mapping_)},
		nodes{knots | std::views::transform(mapping) | std::ranges::to<std::vector>()},
		weights{barycentric_weights(nodes)}
	{
		if (knots.size() < 2) {
			throw std::invalid_argument(std::format("At least two interpolation knots must be supplied, got {}", knots.size()));
		}
		if (!(std::ranges::is_sorted(knots) && knots.front() < knots.back())) {
			throw std::invalid_argument("Interpolation knots must be in increasing order");
		}
		if (!collections::is_unique(knots)) {
			throw std::invalid_argument("Interpolation knots must be unique");
		}
		if (!std::ranges::all_of(nodes, [](const double t) { return std::isfinite(t); })) {
			throw std::invalid_argument("The interpolation mapping must be finite at the knots");
		}
		if (!(std::ranges::is_sorted(nodes) || std::ranges::is_sorted(nodes, std::ranges::greater{}))) {
			throw std::invalid_argument("The interpolation mapping must be monotonic at the knots");
		}
		if (!collections::is_unique(nodes)) {
			throw std::invalid_argument("The interpolation mapping must be injective at the knots");
		}

		kernels = construct_kernels(nodes, weights);
	}

	/// Returns the number of kernels
	std::size_t size() const noexcept { return kernels.size(); }

	/// Returns the kernels that are non-zero at \p x, as a range
	auto nonzero_at(const double x) const {
		using diff_t = typename decltype(kernels)::difference_type;

		const auto [start_idx, end_idx] = [&]() -> std::pair<diff_t, diff_t> {
			const bool in_range = knots.front() <= x && x <= knots.back();

			// Outside the knot range, there are no non-zero kernels
			if (!in_range) {
				return {0, 0};
			}

			const auto node = std::ranges::find(nodes, mapping(x));

			// At a node only that node is non-zero
			if (node != nodes.end()) {
				const diff_t node_idx = collections::iterator_to_index(nodes, node);
				return {node_idx, node_idx + 1};
			}

			// Otherwise, all nodes contribute
			return {0, std::ssize(kernels)};
		}();

		return kernels | shims::views::slice(start_idx, end_idx);
	}

	/// Returns the values of the kernels that are non-zero at \p x (via `nonzero_at`) as a range of (index, value) pairs
	auto values_at(const double x) const {
		const double t = mapping(x);
		const double weight_sum = std::ranges::fold_left(
			std::views::zip_transform([t](const double node, const double weight) { return weight / (t - node); }, nodes, weights), 0.0, std::plus{}
		);
		return nonzero_at(x) | std::views::transform([t, weight_sum](const kernel &k) { return std::pair{k.index, k(t, weight_sum)}; });
	}

private:
	/// Returns the barycentric weights of \p nodes_, which are given by the formula `w_k = 1 / prod_{j != k} (t_k - t_j)`
	static std::vector<double> barycentric_weights(const std::vector<double> &nodes_) {
		return nodes_ | std::views::transform([&nodes_](const double node) {
			return 1.0 / std::ranges::fold_left(
				nodes_
					| std::views::filter([node](const double other) { return other != node; })
					| std::views::transform([node](const double other) { return node - other; }),
				1.0, std::multiplies{}
			);
		}) | std::ranges::to<std::vector>();
	}

	/// Returns the kernels of \p nodes_ with the barycentric weights \p weights_, one for each node
	static std::vector<kernel> construct_kernels(const std::vector<double> &nodes_, const std::vector<double> &weights_) {
		return std::views::zip_transform(
			[](const kernel_index index, const double node, const double weight) { return kernel{index, node, weight}; },
			std::views::iota(kernel_index{0}), nodes_, weights_
		) | std::ranges::to<std::vector>();
	}

	/// The knots, in increasing order
	std::vector<double> knots;
	/// The mapping of the momentum fraction to the variable in which the kernels are polynomials
	mapping_t mapping;
	/// The mapping of the knots
	std::vector<double> nodes;
	/// The barycentric weights of the nodes
	std::vector<double> weights;
	/// The kernels, one for each knot
	std::vector<kernel> kernels;
};

}

#endif
