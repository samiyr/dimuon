#ifndef DM_UTILITY_NAMES_H
#define DM_UTILITY_NAMES_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <format>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace dm {

/// Utilities for providing support for mapping between enumerated values and their string representations 
namespace names {
	/// Type alias for a mapping pair between the enumerated value and the name
	template <typename value_t>
	using entry = std::pair<value_t, std::string_view>;

	/// Essentially a simple bidirectional map - a type alias for a list of mapping pairs
	template <typename value_t, std::size_t size>
	using table = std::array<entry<value_t>, size>;

	/// @brief Returns the name of \p value in \p names
	/// @throws std::out_of_range if \p names has no entry for \p value
	template <typename value_t, std::size_t size>
	constexpr std::string_view name_of(const table<value_t, size> &names, const value_t value) {
		const auto found = std::ranges::find(names, value, &entry<value_t>::first);
		if (found == names.end()) { throw std::out_of_range("The value has no name in the table"); }
		return found->second;
	}

	/// Returns the value by name \p name in \p names, or `std::nullopt` if no value is found by \p name
	template <typename value_t, std::size_t size>
	constexpr std::optional<value_t> value_of(const table<value_t, size> &names, const std::string_view name) noexcept {
		const auto found = std::ranges::find(names, name, &entry<value_t>::second);
		return found != names.end() ? std::optional{found->first} : std::nullopt;
	}

	/// @brief Add formatter specialization to the `std` namespace so that calling 
	/// `std::format` on an enumerated value prints its name as given by `names::name_of`. 
	/// @details Usually adding anything to the `std` namespace is undefined behaviour.
	/// In this instance, however, it is allowed 
	/// (because the specialization depends on a program-defined type)
	template <const auto &names>
	struct formatter : std::formatter<std::string_view> {
		/// Reconstruct `value_t` used in the `names` namespace (i.e. the value type of the first item in the mapping pair)
		using value_t = typename std::remove_cvref_t<decltype(names)>::value_type::first_type;

		auto format(const value_t value, auto &context) const {
			return std::formatter<std::string_view>::format(name_of(names, value), context);
		}
	};
}

}

#endif
