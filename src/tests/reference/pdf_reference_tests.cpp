#include <algorithm>
#include <format>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include <LHAPDF/LHAPDF.h> // IWYU pragma: keep (umbrella header of LHAPDF)

#include <dimuon/common/flavor.hpp>
#include <dimuon/pdf/interface/lha_interface.hpp>

#include "support/assertions.hpp"
#include "support/reference_values.hpp"

using namespace dm;
using namespace dm::testing;

TEST(pdf_reference, lhapdf_values) {
	// The values x f(x, Q^2) of the central member of the PDF set, which tests the flavor numbering and the use of LHAPDF
	const std::string set = "EPPS21nlo_CT18Anlo_Fe56";
	lha_interface<>::disable_verbosity();
	if (!std::ranges::contains(LHAPDF::availablePDFSets(), set)) { GTEST_SKIP() << "The PDF set " << set << " is not installed"; }
	const std::string file = std::format("pdf/{}.tsv", set);
	const std::optional<std::vector<reference_point>> reference = read_reference_points(file);
	if (!reference) { GTEST_SKIP() << missing_reference_message(file); }
	const lha_interface pdf{set};

	for (const auto &[inputs, expected] : *reference) {
		ASSERT_EQ(inputs.size(), 3) << file << ": a row is not \"x Q2 flavor value\"";
		const double x = inputs[0], Q2 = inputs[1];
		const auto parton = std::ranges::find(flavor::all_flavors, inputs[2], [](const flavor::type f) { return static_cast<double>(std::to_underlying(f)); });
		ASSERT_NE(parton, flavor::all_flavors.end()) << file << ": " << inputs[2] << " is not a flavor";
		EXPECT_REL_NEAR(pdf.xfxQ2(*parton, x, Q2), expected, 1e-14) << file << ", flavor = " << inputs[2] << ", x = " << x << ", Q2 = " << Q2;
	}
}
