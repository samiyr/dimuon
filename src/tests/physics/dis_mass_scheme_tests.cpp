#include <utility>

#include <gtest/gtest.h>

#include <dimuon/common/mass_scheme.hpp>
#include <dimuon/common/perturbative_quantity.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/common/structure_function.hpp>
#include <dimuon/common/trf_kinematics.hpp>

#include "support/assertions.hpp"
#include "support/dis_setup.hpp"
#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum scattering_process::process_type;
	using enum structure_function_t;

	const double x = 0.1;
	const double Q2 = 10.0;
	const double s = 300.0;

	auto make_zero_mass_dis(const scattering_process::process_type type) { return make_dis(type, mass_scheme_t::zero_mass, mass_scheme_t::zero_mass, 0.0, toy_pdf_interface()); }

	trf_kinematics make_kinematics(const auto &dis) { return trf_kinematics::Q2_s(x, Q2, s, dis.process.target.mass, dis.process.projectile.mass); }

	/// Expects DIS with a nearly massless charm in the scheme \p massive to agree with the zero-mass scheme, at LO and NLO
	void expect_massless_limit(const mass_scheme_t massive) {
		for (const auto type : {neutrino, antineutrino}) {
			const auto zero_mass = make_zero_mass_dis(type);
			const auto massive_dis = make_dis(type, mass_scheme_t::slow_rescaling, massive, 1e-3, toy_pdf_interface());
			const auto compare = [&]<structure_function_t sf>() {
				const perturbative_quantity a = massive_dis.structure_function<sf>(make_kinematics(zero_mass)).value;
				const perturbative_quantity b = zero_mass.structure_function<sf>(make_kinematics(zero_mass)).value;
				EXPECT_REL_NEAR(a.lo(), b.lo(), 1e-5) << "process " << std::to_underlying(type) << ", F" << std::to_underlying(sf) + 1;
				EXPECT_REL_NEAR(a.nlo(), b.nlo(), 1e-5) << "process " << std::to_underlying(type) << ", F" << std::to_underlying(sf) + 1;
			};
			compare.template operator()<f1>();
			compare.template operator()<f2>();
			compare.template operator()<f3>();
		}
	}
}

TEST(dis_mass_scheme, slow_rescaling_limit) {
	expect_massless_limit(mass_scheme_t::slow_rescaling);
}

TEST(dis_mass_scheme, sacot_limit) {
	expect_massless_limit(mass_scheme_t::sacot_subtracted);
}
