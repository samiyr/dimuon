#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/constants.hpp>
#include <dimuon/common/flavor.hpp>
#include <dimuon/common/parton.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/coupling/electromagnetic_coupling.hpp>

using namespace dm;

namespace {
	/// All flavors except the gluon
	std::vector<flavor::type> quark_flavors() {
		return flavor::all_flavors | std::views::filter([](const flavor::type f) { return !flavor::is_gluon(f); }) | std::ranges::to<std::vector>();
	}

	/// The flavors of \p partons, sorted
	std::vector<flavor::type> sorted_flavors(const std::vector<parton> &partons) {
		std::vector<flavor::type> flavors = partons | std::views::transform(&parton::flavor) | std::ranges::to<std::vector>();
		std::ranges::sort(flavors);
		return flavors;
	}
}

TEST(flavor, charges) {
	// In units of e/3
	for (const flavor::type up_type : {flavor::u, flavor::c, flavor::t}) { 
		EXPECT_EQ(flavor::charge(up_type), 2);
		EXPECT_EQ(flavor::charge(flavor::conjugate_flavor(up_type)), -2);
	}
	for (const flavor::type down_type : {flavor::d, flavor::s, flavor::b}) { 
		EXPECT_EQ(flavor::charge(down_type), -1); 
		EXPECT_EQ(flavor::charge(flavor::conjugate_flavor(down_type)), 1);
	}
	EXPECT_EQ(flavor::charge(flavor::g), 0);
}

TEST(flavor, conjugation) {
	for (const flavor::type f : quark_flavors()) {
		const flavor::type conjugate = flavor::conjugate_flavor(f);

		EXPECT_EQ(flavor::conjugate_flavor(conjugate), f);
		EXPECT_EQ(flavor::is_quark_flavor(f), flavor::is_quark_antiflavor(conjugate));
		EXPECT_EQ(flavor::charge(conjugate), -flavor::charge(f));
		EXPECT_EQ(flavor::abs(conjugate), flavor::abs(f));
		EXPECT_EQ(flavor::is_upper_flavor(conjugate), flavor::is_upper_flavor(f));
	}
	EXPECT_EQ(flavor::conjugate_flavor(flavor::g), flavor::g);
}

TEST(flavor, quark_types) {
	for (const flavor::type f : quark_flavors()) {
		EXPECT_NE(flavor::is_upper_flavor(f), flavor::is_lower_flavor(f)) << "flavor " << std::to_underlying(f);
		EXPECT_EQ(flavor::is_upper_flavor(f), flavor::charge(flavor::abs(f)) == 2) << "flavor " << std::to_underlying(f);
		EXPECT_EQ(flavor::is_lower_flavor(f), flavor::charge(flavor::abs(f)) == -1) << "flavor " << std::to_underlying(f);
	}
	EXPECT_FALSE(flavor::is_upper_flavor(flavor::g));
	EXPECT_FALSE(flavor::is_lower_flavor(flavor::g));
}

TEST(parton_configuration, construction) {
	std::array<double, constants::total_flavors> masses{};

	// The antiquark and the gluon should be ignored
	const std::vector<flavor::type> allowed{flavor::d, flavor::u, flavor::s, flavor::c, flavor::b, flavor::g, flavor::sbar};
	const parton_configuration configuration{allowed, masses};

	using enum flavor::type;
	EXPECT_EQ(sorted_flavors(configuration.allowed_quarks), (std::vector{d, u, s, c, b}));
	EXPECT_EQ(sorted_flavors(configuration.allowed_antiquarks), (std::vector{bbar, cbar, sbar, ubar, dbar}));
	EXPECT_EQ(sorted_flavors(configuration.allowed_lower_quarks), (std::vector{d, s, b}));
	EXPECT_EQ(sorted_flavors(configuration.allowed_upper_quarks), (std::vector{u, c}));
	EXPECT_EQ(sorted_flavors(configuration.allowed_lower_antiquarks), (std::vector{bbar, sbar, dbar}));
	EXPECT_EQ(sorted_flavors(configuration.allowed_upper_antiquarks), (std::vector{cbar, ubar}));
}

TEST(parton_configuration, masses) {
	const double charm_mass = 1.3;
	const double bottom_mass = 4.75;
	std::array<double, constants::total_flavors> masses{};
	masses[flavor::to_index(flavor::c)] = charm_mass;
	masses[flavor::to_index(flavor::b)] = bottom_mass;

	const std::vector<flavor::type> allowed{flavor::d, flavor::u, flavor::s, flavor::c, flavor::b, flavor::g, flavor::sbar};
	const parton_configuration configuration{allowed, masses};

	using enum flavor::type;
	// An antiquark has the mass of its quark
	for (const parton &antiquark : configuration.allowed_antiquarks) {
		EXPECT_EQ(antiquark.mass, masses[flavor::to_index(flavor::conjugate_flavor(antiquark.flavor))]);
	}

	EXPECT_EQ(configuration.mass(c), charm_mass);
	EXPECT_EQ(configuration.mass(cbar), charm_mass);
	EXPECT_EQ(configuration.mass(b), bottom_mass);
	EXPECT_EQ(configuration.mass(u), 0.0);
	EXPECT_EQ(configuration.mass(g), 0.0);

	// The top quark was not included in the list of allowed quarks
	EXPECT_THROW(static_cast<void>(configuration.mass(t)), std::out_of_range);
}

TEST(ckm_matrix, unitarity) {
	// The CKM matrix is unitary (the actual fit values are unitary within the fit tolerance)
	const ckm_matrix V = ckm_matrix::global_fit();

	for (const flavor::type up_type : {flavor::u, flavor::c, flavor::t}) {
		EXPECT_NEAR(V(up_type, flavor::d) + V(up_type, flavor::s) + V(up_type, flavor::b), 1.0, 1e-4) << "row of flavor " << std::to_underlying(up_type);
	}
	for (const flavor::type down_type : {flavor::d, flavor::s, flavor::b}) {
		EXPECT_NEAR(V(flavor::u, down_type) + V(flavor::c, down_type) + V(flavor::t, down_type), 1.0, 1e-4) << "column of flavor " << std::to_underlying(down_type);
	}
}

TEST(ckm_matrix, upper_lower_quark_coupling) {
	const ckm_matrix V = ckm_matrix::global_fit();

	for (const flavor::type a : quark_flavors()) {
		for (const flavor::type b : quark_flavors()) {
			if (flavor::is_upper_flavor(a) != flavor::is_upper_flavor(b)) {
				EXPECT_GT(V(a, b), 0.0) << std::to_underlying(a) << " " << std::to_underlying(b);
				// The same element for either order of the quarks and for their antiquarks
				EXPECT_EQ(V(a, b), V(b, a));
				EXPECT_EQ(V(a, b), V(flavor::conjugate_flavor(a), b));
				EXPECT_EQ(V(a, b), V(a, flavor::conjugate_flavor(b)));
			} else {
				EXPECT_EQ(V(a, b), 0.0) << std::to_underlying(a) << " " << std::to_underlying(b);
			}
		}
	}
}

TEST(ckm_matrix, construction) {
	const double V_ud = 1.0, V_us = 2.0, V_ub = 3.0;
	const double V_cd = 4.0, V_cs = 5.0, V_cb = 6.0;
	const double V_td = 7.0, V_ts = 8.0, V_tb = 9.0;
	const ckm_matrix V{{V_ud, V_us, V_ub, V_cd, V_cs, V_cb, V_td, V_ts, V_tb}};
	using enum flavor::type;
	EXPECT_EQ(V(u, d), V_ud);
	EXPECT_EQ(V(u, s), V_us);
	EXPECT_EQ(V(u, b), V_ub);
	EXPECT_EQ(V(c, d), V_cd);
	EXPECT_EQ(V(c, s), V_cs);
	EXPECT_EQ(V(c, b), V_cb);
	EXPECT_EQ(V(t, d), V_td);
	EXPECT_EQ(V(t, s), V_ts);
	EXPECT_EQ(V(t, b), V_tb);
}

TEST(electromagnetic_coupling, construction) {
	const electromagnetic_coupling coupling;
	const double down_coupling = coupling(flavor::d, flavor::d);
	
	for (const flavor::type a : quark_flavors()) {
		for (const flavor::type b : quark_flavors()) {
			if (flavor::abs(a) == flavor::abs(b)) {
				const double squared_charge_ratio = static_cast<double>(math::power<2>(flavor::charge(a)));
				EXPECT_DOUBLE_EQ(coupling(a, b), squared_charge_ratio * down_coupling) << std::to_underlying(a) << " " << std::to_underlying(b);
			} else {
				EXPECT_EQ(coupling(a, b), 0.0) << std::to_underlying(a) << " " << std::to_underlying(b);
			}
		}
	}
	EXPECT_EQ(coupling(flavor::g, flavor::g), 0.0);
}
