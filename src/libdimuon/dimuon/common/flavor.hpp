#ifndef DM_COMMON_FLAVOR_H
#define DM_COMMON_FLAVOR_H

#include <array>
#include <concepts>
#include <cstdlib>
#include <format>
#include <type_traits>
#include <utility>
#include <vector>

#include "dimuon/utility/names.hpp"

namespace dm {

/// Provides access and convenience methods for handling quark flavors (+ gluon) using the MC particle identification scheme.
/// An exception is made for the gluon, which has the ID 0 instead of 21, since this avoids exceptions for the gluon and is supported
/// by LHAPDF.
namespace flavor {
	/// Parton identities
	enum class type : int {
		u = 2,
		d = 1,
		c = 4,
		s = 3,
		t = 6,
		b = 5,
		g = 0,
		ubar = -2,
		dbar = -1,
		cbar = -4,
		sbar = -3,
		tbar = -6,
		bbar = -5
	};
	/// Concrete underlying type of `flavor::type` (i.e. `int`)
	using underlying_flavor_t = std::underlying_type_t<type>;
	/// Shorthand for `std::vector<type>`
	using flavor_vector = std::vector<type>;

	/// The flavors as constants of this namespace (e.g. `flavor::u`)
	using enum type;

	/// All flavors, in the order of their indices (see `to_index`)
	inline constexpr std::array all_flavors{tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t};

	/// The symbols of the flavors (e.g. "ubar"), see `names`
	inline constexpr auto flavor_symbols = std::to_array<names::entry<type>>({
		{u, "u"},
		{d, "d"},
		{c, "c"},
		{s, "s"},
		{t, "t"},
		{b, "b"},
		{g, "g"},
		{ubar, "ubar"},
		{dbar, "dbar"},
		{cbar, "cbar"},
		{sbar, "sbar"},
		{tbar, "tbar"},
		{bbar, "bbar"}
	});

	/// Applies an arbitrary operation to the underlying value of a flavor and casts the result back to `flavor::type`
	template <typename operation_t>
	requires requires(operation_t operation, const underlying_flavor_t f) {
		{ operation(f) } -> std::convertible_to<underlying_flavor_t>;
	}
	constexpr type apply_operation(const type flavor, operation_t &&operation) noexcept(noexcept(std::forward<operation_t>(operation)(std::to_underlying(flavor)))) {
		return static_cast<type>(std::forward<operation_t>(operation)(std::to_underlying(flavor)));
	}

	/// Converts the given flavor to a corresponding flavor vector index.
	/// This simply adds 6 to the flavor number (so that tbar = -6 becomes index = 0)
	/// and casts the result to std::size_t.
	constexpr std::size_t to_index(const type flavor) noexcept {
		return static_cast<std::size_t>(std::to_underlying(flavor) + 6);
	}

	/// Converts the given flavor to its conjugate flavor, i.e. turn a quark flavor to an antiquark flavor and vice versa.
	/// Does nothing to a gluon.
	constexpr type conjugate_flavor(const type flavor) noexcept {
		return apply_operation(flavor, [](const underlying_flavor_t f) { return -f; });
	}

	constexpr type abs(const type flavor) noexcept {
		return apply_operation(flavor, [](const underlying_flavor_t f) { return std::abs(f); });
	}

	/// Checks whether a given flavor is an upper flavor (i.e. u, c, t). A gluon is not an upper flavor.
	constexpr bool is_upper_flavor(const type flavor) noexcept {
		return std::to_underlying(flavor) != 0 && std::to_underlying(flavor) % 2 == 0;
	}
	/// Checks whether a given flavor is a lower flavor (i.e. d, s, b). A gluon is not a lower flavor.
	constexpr bool is_lower_flavor(const type flavor) noexcept {
		return std::abs(std::to_underlying(flavor)) % 2 == 1;
	}
	/// Checks whether a given flavor is a quark flavor. A gluon is not a quark flavor.
	constexpr bool is_quark_flavor(const type flavor) noexcept {
		return std::to_underlying(flavor) > 0;
	}
	/// Checks whether a given flavor is a quark antiflavor. A gluon is not an antiflavor.
	constexpr bool is_quark_antiflavor(const type flavor) noexcept {
		return std::to_underlying(flavor) < 0;
	}
	/// Checks whether the given flavor is the gluon
	constexpr bool is_gluon(const type flavor) noexcept {
		return flavor == flavor::g;
	}

	/// Gives the electric charge of the given parton flavor, in units of e/3.
	/// For example, charge of u = 2, charge of d = -1.
	/// Returns 0 for the gluon.
	constexpr int charge(const type flavor) noexcept {
		if (is_gluon(flavor)) { return 0; }
		
		const int multiplier = is_quark_antiflavor(flavor) ? -1 : 1;

		if (is_upper_flavor(flavor)) { 
			return multiplier * 2; 
		} else if (is_lower_flavor(flavor)) { 
			return -multiplier * 1; 
		}

		return 0;
	}
};

}

/// Formats a flavor as its symbol (see `dm::flavor::flavor_symbols`)
template <> struct std::formatter<dm::flavor::type> : dm::names::formatter<dm::flavor::flavor_symbols> {};

#endif