#ifndef DM_COMMON_COEFFICIENT_IMPLEMENTATION_H
#define DM_COMMON_COEFFICIENT_IMPLEMENTATION_H

#include <array>
#include <format>

#include "dimuon/utility/names.hpp"

namespace dm {

/// Type of coefficient implementation
enum class coefficient_implementation_t { 
	/// The "handwritten" implementation, i.e.
	/// coefficients that have been converted to code
	/// by hand
	hw,
	/// The "code-generation" implementation, i.e.
	/// coefficients that have been converted to code
	/// with Mathematica
	gc 
};

/// Mapping of coefficient implementations to their corresponding (argument) names, see `names`
inline constexpr auto coefficient_implementation_names = std::to_array<names::entry<coefficient_implementation_t>>({
	{coefficient_implementation_t::hw, "hw"},
	{coefficient_implementation_t::gc, "gc"}
});

}

/// Formats a coefficient implementation as its name (see `dm::coefficient_implementation_names`)
template <> struct std::formatter<dm::coefficient_implementation_t> : dm::names::formatter<dm::coefficient_implementation_names> {};

#endif