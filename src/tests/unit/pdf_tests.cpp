#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include <gtest/gtest.h>
#include <LHAPDF/LHAPDF.h> // IWYU pragma: keep (umbrella header of LHAPDF)

#include <dimuon/common/constants.hpp>
#include <dimuon/common/flavor.hpp>
#include <dimuon/pdf/extrapolator/freeze_extrapolator.hpp>
#include <dimuon/pdf/interface/functional_form_interface.hpp>
#include <dimuon/pdf/interface/lha_interface.hpp>

#include "support/assertions.hpp"

using namespace dm;

namespace {
	using enum flavor::type;

	/// The PDF set of the example in the documentation
	const std::string pdf_set = "EPPS21nlo_CT18Anlo_Fe56";

	/// Points inside the grid of the set
	const std::array<std::pair<double, double>, 4> grid_points{{{1e-4, 2.0}, {0.01, 10.0}, {0.2, 100.0}, {0.7, 1e4}}};

	/// Skips the test if the PDF set is not installed
	class lha_set : public ::testing::Test {
	protected:
		void SetUp() override {
			lha_interface<>::disable_verbosity();
			if (!std::ranges::contains(LHAPDF::availablePDFSets(), pdf_set)) { GTEST_SKIP() << "The PDF set " << pdf_set << " is not installed"; }
		}
	};

	lha_interface<> make_nucleus(const double Z, const double A) {
		lha_interface pdf{pdf_set};
		pdf.enable_isospin(Z, A);
		return pdf;
	}

	void expect_all_flavors_consistent(const auto &pdf) {
		for (const auto &[x, Q2] : grid_points) {
			const std::array<double, constants::total_flavors> values = pdf.xfxQ2(x, Q2);
			for (const flavor::type f : flavor::all_flavors) {
				EXPECT_REL_NEAR(values[flavor::to_index(f)], pdf.xfxQ2(f, x, Q2), 1e-13) << "flavor " << std::to_underlying(f) << " at x = " << x << ", Q2 = " << Q2;
			}
		}
	}
}

TEST(functional_form_interface, all_flavors) {
	const functional_form_interface pdf{[](const flavor::type f, const double x, const double Q2) {
		return static_cast<double>(std::to_underlying(f) + 7) * x * (1.0 - x) * std::log(Q2);
	}};
	expect_all_flavors_consistent(pdf);
}

TEST_F(lha_set, all_flavors) {
	expect_all_flavors_consistent(lha_interface{pdf_set});
	expect_all_flavors_consistent(1.14 * make_nucleus(26.0, 56.0));
}

TEST_F(lha_set, neutron) {
	const lha_interface proton{pdf_set};
	const lha_interface neutron = make_nucleus(0.0, 1.0);

	for (const auto &[x, Q2] : grid_points) {
		EXPECT_REL_NEAR(neutron.xfxQ2(u, x, Q2), proton.xfxQ2(d, x, Q2), 1e-14);
		EXPECT_REL_NEAR(neutron.xfxQ2(d, x, Q2), proton.xfxQ2(u, x, Q2), 1e-14);
		EXPECT_REL_NEAR(neutron.xfxQ2(ubar, x, Q2), proton.xfxQ2(dbar, x, Q2), 1e-14);
		EXPECT_REL_NEAR(neutron.xfxQ2(dbar, x, Q2), proton.xfxQ2(ubar, x, Q2), 1e-14);
	}
}

TEST_F(lha_set, isoscalar) {
	const lha_interface deuteron = make_nucleus(1.0, 2.0);

	for (const auto &[x, Q2] : grid_points) {
		EXPECT_REL_NEAR(deuteron.xfxQ2(u, x, Q2), deuteron.xfxQ2(d, x, Q2), 1e-14);
		EXPECT_REL_NEAR(deuteron.xfxQ2(ubar, x, Q2), deuteron.xfxQ2(dbar, x, Q2), 1e-14);
	}
}

TEST_F(lha_set, nucleus) {
	const double Z = 26.0;
	const double A = 56.0;
	const lha_interface proton{pdf_set};
	const lha_interface iron = make_nucleus(Z, A);
	for (const auto &[x, Q2] : grid_points) {
		EXPECT_REL_NEAR(iron.xfxQ2(u, x, Q2), (Z * proton.xfxQ2(u, x, Q2) + (A - Z) * proton.xfxQ2(d, x, Q2)) / A, 1e-14);
		EXPECT_REL_NEAR(iron.xfxQ2(dbar, x, Q2), (Z * proton.xfxQ2(dbar, x, Q2) + (A - Z) * proton.xfxQ2(ubar, x, Q2)) / A, 1e-14);
		for (const flavor::type f : {s, sbar, c, cbar, b, g}) { EXPECT_EQ(iron.xfxQ2(f, x, Q2), proton.xfxQ2(f, x, Q2)); }
	}
}

TEST_F(lha_set, multiplication) {
	// Multiplied copies leave the original unchanged
	const lha_interface pdf{pdf_set};
	const double x = 0.01;
	const double Q2 = 10.0;
	const double value = pdf.xfxQ2(s, x, Q2);
	EXPECT_REL_NEAR((2.0 * pdf).xfxQ2(s, x, Q2), 2.0 * value, 1e-15);
	EXPECT_REL_NEAR((pdf * 3.0).xfxQ2(x, Q2)[flavor::to_index(s)], 3.0 * value, 1e-15);
	EXPECT_EQ(pdf.xfxQ2(s, x, Q2), value);
}

TEST_F(lha_set, in_place_multiplication) {
	lha_interface pdf{pdf_set};
	const double x = 0.01;
	const double Q2 = 10.0;
	const double value = pdf.xfxQ2(s, x, Q2);
	pdf *= 0.5;
	EXPECT_REL_NEAR(pdf.xfxQ2(s, x, Q2), 0.5 * value, 1e-15);
}

TEST_F(lha_set, zero_extrapolation) {
	const std::unique_ptr<LHAPDF::PDF> grid{LHAPDF::mkPDF(pdf_set, 0)};
	const lha_interface pdf{pdf_set};
	EXPECT_EQ(pdf.xfxQ2(d, 0.1 * grid->xMin(), 10.0), 0.0);
	EXPECT_EQ(pdf.xfxQ2(d, 0.1, 0.5 * grid->q2Min()), 0.0);
	EXPECT_EQ(pdf.xfxQ2(d, 0.1, 2.0 * grid->q2Max()), 0.0);
}

TEST_F(lha_set, freeze_extrapolation) {
	const lha_interface<freeze_extrapolator> frozen{pdf_set};
	const std::unique_ptr<LHAPDF::PDF> grid{LHAPDF::mkPDF(pdf_set, 0)};

	const double x_min = grid->xMin();
	const double Q2_min = grid->q2Min();
	const double Q2_max = grid->q2Max();
	
	const double x = 0.1;
	const double Q2 = 10.0;

	EXPECT_EQ(frozen.xfxQ2(d, 0.1 * x_min, Q2), grid->xfxQ2(1, x_min, Q2));
	EXPECT_EQ(frozen.xfxQ2(d, x, 0.5 * Q2_min), grid->xfxQ2(1, x, Q2_min));
	EXPECT_EQ(frozen.xfxQ2(d, x, 2.0 * Q2_max), grid->xfxQ2(1, x, Q2_max));
	EXPECT_EQ(frozen.xfxQ2(d, 0.1 * x_min, 0.5 * Q2_min), grid->xfxQ2(1, x_min, Q2_min));
}
