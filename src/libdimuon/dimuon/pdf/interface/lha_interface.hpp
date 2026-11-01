#ifndef DM_PDF_INTERFACE_LHA_INTERFACE_H
#define DM_PDF_INTERFACE_LHA_INTERFACE_H

#include <string>
#include <vector>
#include <array>
#include <concepts>
#include <cstddef>
#include <memory>
#include <utility>

#include "LHAPDF/GridPDF.h"

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"

#include "dimuon/pdf/extrapolator/zero_extrapolator.hpp"

namespace dm {

/// @brief LHAPDF wrapper
/// @tparam extrapolator_t Type of the extrapolator class
template <std::derived_from<LHAPDF::Extrapolator> extrapolator_t = zero_extrapolator>
class lha_interface {
	/// Set name
	std::string set_name;
	/// Set member
	std::size_t set_member_number;
public:

	/// Per-flavor multipliers of one for all flavors, the default of `flavor_multipliers`
	static constexpr std::array<double, constants::total_flavors> unit_flavor_multipliers = []() {
		std::array<double, constants::total_flavors> multipliers{};
		multipliers.fill(1.0);
		return multipliers;
	}();

	/// A global multiplier, by which all distribution values are multiplied. Defaults to one.
	double global_multiplier = 1.0;
	/// Per-flavor multipliers, indexed by `flavor::to_index`, by which the values of each flavor are multiplied on top of the global multiplier.
	/// Default to one.
	std::array<double, constants::total_flavors> flavor_multipliers;

	/// @brief Constructs the LHAPDF wrapper
	/// @param[in] set_name_ Set name
	/// @param[in] set_member_number_ Set member (defaults to 0, i.e. the central set)
	/// @param[in] flavor_multipliers_ Per-flavor multipliers, indexed by `flavor::to_index` (defaults to one for all flavors)
	explicit lha_interface(
		std::string set_name_,
		std::size_t set_member_number_ = 0,
		const std::array<double, constants::total_flavors> &flavor_multipliers_ = unit_flavor_multipliers
	) : set_name{std::move(set_name_)}, 
	set_member_number{set_member_number_}, 
	flavor_multipliers{flavor_multipliers_},
	pdf{LHAPDF::mkPDF(set_name, set_member_number)} {
		LHAPDF::Extrapolator *extrapolator = new extrapolator_t();
		// Ownership of extrapolator is passed to LHAPDF::GridPDF, it handles deletion.
		std::static_pointer_cast<LHAPDF::GridPDF>(pdf)->setExtrapolator(extrapolator);
	}

	/// Disables LHAPDF output
	static void disable_verbosity() {
		LHAPDF::Info &config = LHAPDF::getConfig();
		config.set_entry("Verbosity", 0);
	}

	/// Returns the set name
	std::string get_set_name() const { return set_name; }
	/// Returns the set member 
	std::size_t get_set_member() const noexcept { return set_member_number; }

	/// @brief Enables isospin corrections
	/// @param[in] Z Atomic number
	/// @param[in] A Mass number
	void enable_isospin(const double Z, const double A) noexcept {
		isospin_Z = Z;
		isospin_A = A;
		explicit_isospin = true;
	}
	/// Disables isospin corrections
	void disable_isospin() noexcept {
		isospin_Z = 1.0;
		isospin_A = 1.0;
		explicit_isospin = false;
	}

	/// @brief Returns the grid value, taking the multipliers and isospin corrections into account (if enabled)
	/// @param[in] flavor Parton flavor
	/// @param[in] x Momentum fraction
	/// @param[in] Q2 Momentum scale (GeV^2)
	double xfxQ2(const flavor::type flavor, const double x, const double Q2) const {
		return adjusted_value(flavor, [&](const flavor::type f) { return pdf->xfxQ2(std::to_underlying(f), x, Q2); });
	}

	/// @brief Computes the values of all flavors at once, indexed by `flavor::to_index`. 
	/// Gives the same values as calling `xfxQ2` for each flavor separately, 
	/// but is faster, since the interpolation work can be shared among the flavors.
	/// @param[in] x Momentum fraction
	/// @param[in] Q2 Momentum scale (GeV^2)
	std::array<double, constants::total_flavors> xfxQ2(const double x, const double Q2) const {
		pdf->xfxQ2(x, Q2, all_flavor_values);

		std::array<double, constants::total_flavors> values;
		for (const flavor::type flavor : flavor::all_flavors) {
			values[flavor::to_index(flavor)] = adjusted_value(flavor, [&](const flavor::type f) { return all_flavor_values[flavor::to_index(f)]; });
		}
		return values;
	}
	
	/// Returns the `alpha_s` value at \p Q2
	double alpha_s(const double Q2) const { return pdf->alphasQ2(Q2); }

	/// Returns the lower edge of the `Q^2` grid
	double Q2_min() const { return pdf->q2Min(); }
	/// Returns the upper edge of the `Q^2` grid
	double Q2_max() const { return pdf->q2Max(); }
	/// Returns the mass of a quark of flavor \p flavor
	double quark_mass(const flavor::type flavor) const { return pdf->quarkMass(std::to_underlying(flavor)); }
	/// Returns the output dimension (which is 1, since this class returns scalars)
	std::size_t dimension() const noexcept { return 1; }

private:
	std::shared_ptr<LHAPDF::PDF> pdf;

	/// Buffer for the values given by LHAPDF for all flavors, reused to avoid allocating at each all-flavor call of `xfxQ2`
	mutable std::vector<double> all_flavor_values;

	double isospin_Z = 1.0;
	double isospin_A = 1.0;
	bool explicit_isospin = false;

	/// Applies the isospin transformation (if enabled), the global multiplier, and the multiplier of \p flavor to the PDF value of \p flavor.
	/// \p raw_value gives the value of a flavor from the underlying PDF instance.
	double adjusted_value(const flavor::type flavor, const auto &raw_value) const {
		double value = raw_value(flavor);

		if (explicit_isospin) {
			const double Z_factor = isospin_Z / isospin_A;
			const double N_factor = (isospin_A - isospin_Z) / isospin_A;

			bool u_or_d = true;
			double conjugate_value;

			switch (flavor) {
			case flavor::u:
				conjugate_value = raw_value(flavor::d);
				break;

			case flavor::ubar:
				conjugate_value = raw_value(flavor::dbar);
				break;

			case flavor::d:
				conjugate_value = raw_value(flavor::u);
				break;

			case flavor::dbar:
				conjugate_value = raw_value(flavor::ubar);
				break;

			default:
				u_or_d = false;
				break;
			}

			if (u_or_d) {
				value = Z_factor * value + N_factor * conjugate_value;
			}
		}

		value *= global_multiplier * flavor_multipliers[flavor::to_index(flavor)];

		return value;
	}
};

template <std::derived_from<LHAPDF::Extrapolator> extrapolator_t>
lha_interface<extrapolator_t> &operator*=(lha_interface<extrapolator_t> &lhs, const double rhs) noexcept {
	lhs.global_multiplier *= rhs;
	return lhs;
}

template <std::derived_from<LHAPDF::Extrapolator> extrapolator_t> 
lha_interface<extrapolator_t> operator*(const double lhs, const lha_interface<extrapolator_t> &rhs) {
	lha_interface<extrapolator_t> copy{rhs};
	copy *= lhs;
	return copy;
}

template <std::derived_from<LHAPDF::Extrapolator> extrapolator_t> 
lha_interface<extrapolator_t> operator*(const lha_interface<extrapolator_t> &lhs, const double rhs) {
	return rhs * lhs;
}

}

#endif