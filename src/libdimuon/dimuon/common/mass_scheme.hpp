#ifndef DM_COMMON_MASS_SCHEME_H
#define DM_COMMON_MASS_SCHEME_H

#include <array>
#include <format>
#include <utility>

#include "dimuon/utility/names.hpp"

namespace dm {

/// Enumeration of available mass schemes
enum class mass_scheme_t {
	/// No kinematical (convolutions based on Bjorken-x) or dynamical (no mass terms in the coefficients) mass corrections
	zero_mass, 
	/// Only kinematical (convolutions based on the slow-rescaling variable), no dynamical (no mass terms in the coefficients) mass corrections
	slow_rescaling, 
	/// Kinematical (convolutions based on the slow-rescaling variable) and dynamical (mass terms in the coefficients) mass corrections.
	/// Subtraction terms are not included in the coefficients.
	sacot_unsubtracted, 
	/// Kinematical (convolutions based on the slow-rescaling variable) and dynamical (mass terms in the coefficients) mass corrections.
	/// Subtraction terms are included in the coefficients.
	sacot_subtracted,
	/// Only kinematical (convolutions based on the slow-rescaling variable), no dynamical (no mass terms in the coefficients) mass corrections.
	/// The F2 structure function is proportional to only Bjorken-x, not the slow-rescaling variable.
	slow_rescaling_modified, 
	/// Kinematical (convolutions based on the slow-rescaling variable) and dynamical (mass terms in the coefficients) mass corrections.
	/// Subtraction terms are not included in the coefficients.
	/// The F2 structure function is proportional to only Bjorken-x, not the slow-rescaling variable.
	sacot_unsubtracted_modified, 
	/// Kinematical (convolutions based on the slow-rescaling variable) and dynamical (mass terms in the coefficients) mass corrections.
	/// Subtraction terms are included in the coefficients.
	/// The F2 structure function is proportional to only Bjorken-x, not the slow-rescaling variable.
	sacot_subtracted_modified
};

/// Mapping of mass schemes to their corresponding (argument) names, see `names`
inline constexpr auto mass_scheme_names = std::to_array<names::entry<mass_scheme_t>>({
	{mass_scheme_t::zero_mass, "zm"},
	{mass_scheme_t::slow_rescaling, "sr"},
	{mass_scheme_t::sacot_unsubtracted, "gmu"},
	{mass_scheme_t::sacot_subtracted, "gm"},
	{mass_scheme_t::slow_rescaling_modified, "sr*"},
	{mass_scheme_t::sacot_unsubtracted_modified, "gmu*"},
	{mass_scheme_t::sacot_subtracted_modified, "gm*"}
});

/// Returns true if \p scheme is one of the SACOT schemes, whose coefficient functions have mass terms from NLO on.
/// At LO, they are the same as those of the slow-rescaling schemes.
constexpr bool is_sacot_scheme(const mass_scheme_t scheme) noexcept {
	switch (scheme) {
		using enum mass_scheme_t;
		case zero_mass:
		case slow_rescaling:
		case slow_rescaling_modified:
			return false;
		case sacot_unsubtracted:
		case sacot_subtracted:
		case sacot_unsubtracted_modified:
		case sacot_subtracted_modified:
			return true;
	}
	std::unreachable();
}

}

/// Formats a mass scheme as its name (see `dm::mass_scheme_names`)
template <> struct std::formatter<dm::mass_scheme_t> : dm::names::formatter<dm::mass_scheme_names> {};

#endif