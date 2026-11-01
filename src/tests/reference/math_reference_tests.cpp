#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/utility/math.hpp>

#include "support/reference_values.hpp"

using namespace dm;
using namespace dm::testing;

TEST(math_reference, invtanint) {
	const std::optional<std::vector<reference_point>> reference = read_reference_points("math/invtanint.tsv");
	if (!reference) { GTEST_SKIP() << missing_reference_message("math/invtanint.tsv"); }

	for (const auto &[input, expected] : *reference) {
		const double x = input[0];
		EXPECT_NEAR(math::invtanint(x), expected, 2e-15 * std::max(1.0, std::abs(expected))) << "x = " << x;
	}
}
