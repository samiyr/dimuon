#ifndef DM_COMMON_PARTON_CONFIGURATION_H
#define DM_COMMON_PARTON_CONFIGURATION_H

#include <ranges>
#include <unordered_map>
#include <vector>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"
#include "dimuon/common/parton.hpp"

#include "dimuon/utility/collections.hpp"

namespace dm {

/// Structure for storing various allowed parton types
struct parton_configuration {
	/// List of allowed quarks. Does not contain antiquarks or the gluon.
	const std::vector<parton> allowed_quarks;
	/// List of allowed antiquarks, computed from `allowed_quarks`
	const std::vector<parton> allowed_antiquarks;

	/// List of allowed lower quarks (not antiquarks), computed from `allowed_quarks`
	const std::vector<parton> allowed_lower_quarks;
	/// List of allowed upper quarks (not antiquarks), computed from `allowed_quarks`
	const std::vector<parton> allowed_upper_quarks;
	/// List of allowed lower antiquarks, computed from `allowed_antiquarks`
	const std::vector<parton> allowed_lower_antiquarks;
	/// List of allowed upper antiquarks, computed from `allowed_antiquarks`
	const std::vector<parton> allowed_upper_antiquarks;

	/// @brief Constructs with the given list of allowed partons
	/// @param[in] allowed_quarks_ List of allowed quarks. Antiquarks and the gluon will be ignored.
	explicit parton_configuration(collections::range_of_value<parton> auto allowed_quarks_)
	: allowed_quarks{
		std::move(allowed_quarks_) | std::views::filter([](const parton &p) { return flavor::is_quark_flavor(p.flavor); }) | std::ranges::to<std::vector>()
	},
	allowed_antiquarks{parton_configuration::construct_allowed_antiquarks(allowed_quarks)},
	allowed_lower_quarks{parton_configuration::construct_allowed_lower_flavors(allowed_quarks)},
	allowed_upper_quarks{parton_configuration::construct_allowed_upper_flavors(allowed_quarks)},
	allowed_lower_antiquarks{parton_configuration::construct_allowed_lower_flavors(allowed_antiquarks)},
	allowed_upper_antiquarks{parton_configuration::construct_allowed_upper_flavors(allowed_antiquarks)},
	masses{parton_configuration::construct_mass_map(allowed_quarks)}
	{ }

	/// @brief Constructs with the given allowed quark flavors and masses
	/// @param[in] allowed_quark_flavors List of allowed quark flavors. Antiquark flavors and the gluon will be ignored.
	/// @param[in] masses_ An array of quark masses in the order tbar, bbar, ..., gluon, d, u, ..., t. Must contain all 13 entries, even if not all are allowed flavors.
	explicit parton_configuration(std::span<const flavor::type> allowed_quark_flavors, std::span<const double, constants::total_flavors> masses_)
	: parton_configuration{
		allowed_quark_flavors | std::views::transform(
			[masses_](const flavor::type flavor) { 
				return parton { flavor, masses_[flavor::to_index(flavor)] };
			}
		)
	} {}
		
	/// Returns the assigned quark mass for the given flavor
	double mass(const flavor::type flavor) const {
		return masses.at(flavor::abs(flavor));
	}

private:
	/// A dictionary of parton flavors and masses. The keys consist of the flavors in `allowed_quarks`, and of the gluon (with zero mass).
	const std::unordered_map<flavor::type, double> masses;

	/// Constructs the mass dictionary from the given range of allowed quarks
	static std::unordered_map<flavor::type, double> construct_mass_map(const std::ranges::range auto &allowed_quarks) {
		std::unordered_map<flavor::type, double> map;
		for (const parton &p : allowed_quarks) {
			map[p.flavor] = p.mass;
		}
		map[flavor::g] = 0.0;
		return map;
	}

	/// Constructs the vector of allowed antiquarks from the range of allowed quarks by conjugating all flavors
	static std::vector<parton> construct_allowed_antiquarks(const std::ranges::range auto &allowed_quarks) {
		return allowed_quarks | std::views::transform(parton::conjugate_parton) | std::ranges::to<std::vector>();
	}
	/// Constructs the vector of allowed lower quark flavors by filtering out any non-lower-flavored partons using `flavor::is_lower_flavor`
	static std::vector<parton> construct_allowed_lower_flavors(const std::ranges::range auto &allowed_quarks) {
		return allowed_quarks | std::views::filter([](const parton &p) { return flavor::is_lower_flavor(p.flavor); }) | std::ranges::to<std::vector>();
	}
	/// Constructs the vector of allowed upper quark flavors by filtering out any non-upper-flavored partons using `flavor::is_upper_flavor`
	static std::vector<parton> construct_allowed_upper_flavors(const std::ranges::range auto &allowed_quarks) {
		return allowed_quarks | std::views::filter([](const parton &p) { return flavor::is_upper_flavor(p.flavor); }) | std::ranges::to<std::vector>();
	}
};

}

#endif