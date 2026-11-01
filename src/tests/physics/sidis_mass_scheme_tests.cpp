#include <format>
#include <utility>

#include <gtest/gtest.h>

#include <dimuon/common/structure_function.hpp>

#include "support/assertions.hpp"
#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum structure_function_t;

	/// A run with every channel in the scheme \p scheme and nearly massless charm and bottom quarks
	sidis_run make_massless_run(const mass_scheme_t scheme) {
		return sidis_run{.massless_scheme = scheme, .massive_scheme = scheme, .massive_scheme_nnlo = scheme, .charm_mass = 1e-4, .bottom_mass = 1e-4};
	}
}

TEST(sidis_mass_scheme, slow_rescaling_limit) {
	// The zero-mass coefficient functions at the slow-rescaling variable, which becomes x for massless quarks
	const sidis::result zero_mass = compute_sidis(make_massless_run(mass_scheme_t::zero_mass), toy_pdf_interface(), toy_ff_interface());
	const sidis::result slow_rescaling = compute_sidis(make_massless_run(mass_scheme_t::slow_rescaling), toy_pdf_interface(), toy_ff_interface());
	for (const perturbative_order pto : {lo, nlo, nnlo}) {
		EXPECT_REL_NEAR(slow_rescaling.value.correction(pto), zero_mass.value.correction(pto), 1e-7) << "order " << std::to_underlying(pto);
	}
}

TEST(sidis_mass_scheme, sacot_limit) {
	for (const sidis_run &point : {sidis_run{}, sidis_run{.x = 0.3, .y = 0.9, .renormalization = 2.0, .factorization = 0.5, .fragmentation = 3.0}}) {
		sidis_run zero_mass = point;
		zero_mass.order = nlo;
		zero_mass.massless_scheme = zero_mass.massive_scheme = mass_scheme_t::zero_mass;
		zero_mass.charm_mass = zero_mass.bottom_mass = 1e-2;
		zero_mass.method = integration_method::cuba_vegas;
		zero_mass.points = 500'000;
		sidis_run sacot = zero_mass;
		sacot.massive_scheme = mass_scheme_t::sacot_subtracted;

		const auto compare = [&]<structure_function_t sf>() {
			const auto compute = [&](const sidis_run &run) {
				return compute_sidis(run, toy_pdf_interface(), toy_ff_interface(), structure_function::observable<sf>{}).value;
			};
			const perturbative_quantity a = compute(sacot), b = compute(zero_mass);
			EXPECT_REL_NEAR(a.lo_correction(), b.lo_correction(), 1e-6) << std::format("F{} at x = {}", std::to_underlying(sf) + 1, point.x);
			EXPECT_REL_NEAR(a.nlo_correction(), b.nlo_correction(), 1e-2) << std::format("F{} at x = {}", std::to_underlying(sf) + 1, point.x);
		};
		compare.template operator()<f1>();
		compare.template operator()<f2>();
		compare.template operator()<f3>();
	}
}

TEST(sidis_coefficients, implementations) {
	for (const mass_scheme_t scheme : {mass_scheme_t::zero_mass, mass_scheme_t::slow_rescaling, mass_scheme_t::slow_rescaling_modified}) {
		const sidis_run hw{.order = nlo, .massless_scheme = scheme, .massive_scheme = scheme};
		sidis_run gc = hw;
		gc.nlo_implementation = coefficient_implementation_t::gc;
		const sidis::result a = compute_sidis(hw, toy_pdf_interface(), toy_ff_interface());
		const sidis::result b = compute_sidis(gc, toy_pdf_interface(), toy_ff_interface());
		EXPECT_REL_NEAR(a.value.nlo_correction(), b.value.nlo_correction(), 1e-14) << names::name_of(mass_scheme_names, scheme);
	}
}
