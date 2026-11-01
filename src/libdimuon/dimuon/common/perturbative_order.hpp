#ifndef DM_COMMON_PERTURBATIVE_ORDER_H
#define DM_COMMON_PERTURBATIVE_ORDER_H

#include <array>
#include <format>

#include "dimuon/utility/names.hpp"

namespace dm {

/// The perturbative order, from leading order (LO) to next-to-next-to-leading order (NNLO)
enum class perturbative_order {
	lo = 0, nlo = 1, nnlo = 2
};

/// Mapping of perturbative orders to their names, see `names`
inline constexpr auto perturbative_order_names = std::to_array<names::entry<perturbative_order>>({
	{perturbative_order::lo, "LO"},
	{perturbative_order::nlo, "NLO"},
	{perturbative_order::nnlo, "NNLO"}
});

}

/// Formats a perturbative order as its name (see `dm::perturbative_order_names`)
template <> struct std::formatter<dm::perturbative_order> : dm::names::formatter<dm::perturbative_order_names> {};

#endif
