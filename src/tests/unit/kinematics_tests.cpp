#include <cmath>
#include <numbers>
#include <stdexcept>

#include <gtest/gtest.h>

#include <dimuon/common/constants.hpp>
#include <dimuon/common/cross_section_functions.hpp>
#include <dimuon/common/cross_section_observable.hpp>
#include <dimuon/common/perturbative_quantity.hpp>
#include <dimuon/common/scale_dependence.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/common/trf_kinematics.hpp>
#include <dimuon/utility/math.hpp>

#include "support/assertions.hpp"

using namespace dm;

namespace {
	using enum scattering_process::process_type;

	const double proton_mass = constants::particles::proton.mass;

	scattering_process make_process(const scattering_process::process_type type, const particle target = constants::particles::proton) {
		return scattering_process{type, target, constants::particles::neutrino};
	}

	trf_kinematics make_kinematics(const double x, const double y, const double E) {
		return trf_kinematics::y_E_beam(x, y, E, proton_mass, 0.0);
	}
}

TEST(trf_kinematics, parametrizations) {
	// The same kinematics given by (y, E) or by (Q^2, s) or (Q^2, sqrt(s))
	for (const double y : {0.1, 0.9}) {
		for (const double E : {20.0, 300.0}) {
			const double x = 0.3;
			const trf_kinematics from_beam = make_kinematics(x, y, E);
			const trf_kinematics from_s = trf_kinematics::Q2_s(x, from_beam.Q2, from_beam.s, proton_mass, 0.0);
			const trf_kinematics from_sqrt_s = trf_kinematics::Q2_sqrt_s(x, from_beam.Q2, std::sqrt(from_beam.s), proton_mass, 0.0);

			EXPECT_REL_NEAR(from_s.y, y, 1e-14);
			EXPECT_REL_NEAR(from_s.E_beam, E, 1e-14);
			EXPECT_REL_NEAR(from_sqrt_s.y, y, 1e-14);
			EXPECT_REL_NEAR(from_sqrt_s.E_beam, E, 1e-14);
		}
	}
}

TEST(trf_kinematics, virtuality) {
	const double x = 0.3;
	const double y = 0.6;
	const double E = 100.0;
	// Q^2 = x y (s - M^2) for a massless projectile
	const trf_kinematics kinematics = make_kinematics(x, y, E);
	EXPECT_REL_NEAR(kinematics.Q2, x * y * (kinematics.s - proton_mass * proton_mass), 1e-14);
	// Q^2 = 2 M E x y in the TRF
	EXPECT_REL_NEAR(kinematics.Q2, 2.0 * proton_mass * E * x * y, 1e-14);
}

TEST(cross_section_observable, linearity) {
	const trf_kinematics kinematics = make_kinematics(0.2, 0.4, 100.0);
	const auto observable = [&](const double f1, const double f2, const double f3) { return cross_section_observable{}(kinematics, make_process(neutrino), f1, f2, f3); };

	EXPECT_EQ(observable(0.0, 0.0, 0.0), 0.0);

	const std::vector<double> values = math::linear_space(1.0, 10.0, 1.0);
	const double f1 = 1.0;

	for (const double v1 : values) {
		for (const double v2 : values) {
			for (const double v3 : values) {
				EXPECT_REL_NEAR(observable(v1, v2, v3) + observable(f1, v2, v1), observable(f1 + v1, 2.0 * v2, v1 + v3), 1e-12);
				for (const double v4 : math::linear_space(-5.0, 5.0, 1.0)) {
					EXPECT_REL_NEAR(v4 * observable(v1, v2, v3), observable(v4 * v1, v4 * v2, v4 * v3), 1e-12);
				}
			}
		}
	}
}

TEST(cross_section_observable, antineutrino) {
	// F3 has the opposite sign for antineutrinos
	const trf_kinematics kinematics = make_kinematics(0.2, 0.4, 100.0);
	const cross_section_observable observable;
	const double f1 = 1.0;
	const double f2 = 2.0;
	const double f3 = 3.0;
	EXPECT_DOUBLE_EQ(observable(kinematics, make_process(neutrino), f1, f2, f3), observable(kinematics, make_process(antineutrino), f1, f2, -f3));
}

TEST(cross_section_observable, electromagnetic) {
	// F3 doesn't contribute to parity-conserving processes
	const trf_kinematics kinematics = make_kinematics(0.2, 0.4, 100.0);
	const cross_section_observable observable;
	const double f1 = 1.0;
	const double f2 = 2.0;
	EXPECT_DOUBLE_EQ(observable(kinematics, make_process(electromagnetic), f1, f2, 0.0), observable(kinematics, make_process(electromagnetic), f1, f2, 100.0));
}

TEST(cross_section_functions, modified_prefactor) {
	const double E = 200.0;
	const trf_kinematics kinematics = make_kinematics(0.2, 0.6, E);
	const double conversion = math::power<2>(constants::fermi_coupling) * proton_mass * E / std::numbers::pi;
	EXPECT_REL_NEAR(cross_section_functions::prefactor(kinematics, true), 100.0 * cross_section_functions::prefactor(kinematics, false) / conversion, 1e-12);
}

TEST(perturbative_quantity, orders) {
	const double lo = 1.0;
	const double nlo = -2.0;
	const double nnlo = 3.0;
	const double alpha_s = 0.3;

	const double a = alpha_s / (2.0 * std::numbers::pi);
	const perturbative_quantity q{alpha_s, lo, nlo, nnlo};

	EXPECT_DOUBLE_EQ(q.lo(), lo);
	EXPECT_DOUBLE_EQ(q.nlo(), lo + a * nlo);
	EXPECT_DOUBLE_EQ(q.nnlo(), lo + a * nlo + math::power<2>(a) * nnlo);
	EXPECT_DOUBLE_EQ(q.nlo_correction_with_alpha_s(), a * nlo);
	EXPECT_DOUBLE_EQ(q.nnlo_correction_with_alpha_s(), math::power<2>(a) * nnlo);
}

TEST(perturbative_quantity, single_order) {
	const double nlo = 4.0;
	const perturbative_quantity q{perturbative_order::nlo, nlo};
	EXPECT_EQ(q.lo_correction(), 0.0);
	EXPECT_EQ(q.nlo_correction(), nlo);
	EXPECT_EQ(q.nnlo_correction(), 0.0);
}

TEST(perturbative_quantity, arithmetic) {
	const double lo1 = 1.0;
	const double nlo1 = -2.0;
	const double nnlo1 = 3.0;
	const double lo2 = -4.0;
	const double nlo2 = 5.0;
	const double nnlo2 = -6.0;

	const double alpha_s = 0.3;

	const perturbative_quantity p{alpha_s, lo1, nlo1, nnlo1};
	const perturbative_quantity q{alpha_s, lo2, nlo2, nnlo2};

	const auto expect = [](const perturbative_quantity &actual, const double lo, const double nlo, const double nnlo) {
		EXPECT_DOUBLE_EQ(actual.lo_correction(), lo);
		EXPECT_DOUBLE_EQ(actual.nlo_correction(), nlo);
		EXPECT_DOUBLE_EQ(actual.nnlo_correction(), nnlo);
	};

	expect(p + q, lo1 + lo2, nlo1 + nlo2, nnlo1 + nnlo2);
	expect(p - q, lo1 - lo2, nlo1 - nlo2, nnlo1 - nnlo2);
	expect(p * q, lo1 * lo2, nlo1 * nlo2, nnlo1 * nnlo2);
	expect(p / q, lo1 / lo2, nlo1 / nlo2, nnlo1 / nnlo2);
	expect(2.0 * p, 2.0 * lo1, 2.0 * nlo1, 2.0 * nnlo1);
	expect(p * 2.0, 2.0 * lo1, 2.0 * nlo1, 2.0 * nnlo1);
}

TEST(perturbative_quantity, mismatched_alpha_s) {
	const perturbative_quantity p{0.2, 1.0, 2.0, 3.0};
	const perturbative_quantity q{0.3, 1.0, 2.0, 3.0};
	EXPECT_THROW(p + q, std::runtime_error);
	EXPECT_THROW(p - q, std::runtime_error);
	EXPECT_THROW(p * q, std::runtime_error);
	EXPECT_THROW(p / q, std::runtime_error);
}

TEST(scale_dependence, multiplier) {
	const scale_dependence::function scale{[](const trf_kinematics &kinematics) { return kinematics.Q2 + 1.69; }};
	const trf_kinematics kinematics = trf_kinematics::Q2_s(0.1, 10.0, 1000.0, proton_mass, 0.0);
	EXPECT_DOUBLE_EQ((4.0 * scale)(kinematics), 4.0 * scale(kinematics));
	EXPECT_DOUBLE_EQ((scale * 4.0)(kinematics), 4.0 * scale(kinematics));
	EXPECT_DOUBLE_EQ((0.25 * scale)(kinematics), 0.25 * scale(kinematics));
}

TEST(scale_dependence, lower_bound) {
	const double bound = 1.69;
	const scale_dependence::function scale{[](const trf_kinematics &kinematics) { return kinematics.Q2 + 1.69; }, bound};
	const trf_kinematics kinematics = trf_kinematics::Q2_s(0.1, 0.5, 1000.0, proton_mass, 0.0);
	EXPECT_DOUBLE_EQ((0.25 * scale)(kinematics), bound);
}
