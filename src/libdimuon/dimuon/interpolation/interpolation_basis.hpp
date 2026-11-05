#ifndef DM_INTERPOLATION_INTERPOLATION_BASIS_H
#define DM_INTERPOLATION_INTERPOLATION_BASIS_H

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dimuon/utility/collections.hpp"

namespace dm {

/// Index of an interpolation kernel
using kernel_index = std::size_t;

template <typename T>
concept interpolation_kernel_mapping = requires(T mapping, const double x) {
	{ mapping(x) } -> std::convertible_to<double>;
};

template <typename T>
concept interpolation_basis = interpolation_kernel_mapping<typename T::mapping_type>
	&& requires(T basis, const double x, const std::vector<double> knots, const typename T::mapping_type &mapping, const kernel_index index_offset) {
	{ basis.size() } -> std::convertible_to<std::size_t>;
	{ basis.in_range(x) } -> std::same_as<bool>;
	{ basis.nonzero_at(x) } -> collections::range_of_value<typename T::kernel>;
	{ basis.values_at(x) } -> collections::range_convertible_to<std::pair<kernel_index, double>>;
	{ T{knots, mapping, index_offset} };
	{ T::degree } -> std::convertible_to<std::size_t>;
	{ T::maximum_kernels_per_point } -> std::convertible_to<std::size_t>;
};

/// @brief Checks the knots of an interpolation basis and their nodes (the knots mapped by the interpolation mapping), and throws `std::invalid_argument` if they are invalid
/// @param[in] knots The knots, which must be at least two, in increasing order and unique
/// @param[in] nodes The nodes, which must be finite, monotonic and unique, i.e. the mapping must be finite and strictly monotonic at the knots
inline void check_interpolation_knots(const std::vector<double> &knots, const std::vector<double> &nodes) {
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
}

}

#endif