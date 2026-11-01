#ifndef DM_TESTS_SUPPORT_ASSERTIONS_H
#define DM_TESTS_SUPPORT_ASSERTIONS_H

#include <algorithm>
#include <cmath>
#include <format>

#include <gtest/gtest.h>

namespace dm::testing {
	/// Predicate formatter of `EXPECT_REL_NEAR`: \p a and \p b are near if `|a - b| <= tolerance * max(|a|, |b|)`.
	/// Two zeros are near, and a NaN is near nothing.
	inline ::testing::AssertionResult relative_near(
		const char *a_expression, const char *b_expression, const char *tolerance_expression,
		const double a, const double b, const double tolerance
	) {
		const double scale = std::max(std::abs(a), std::abs(b));
		if (std::abs(a - b) <= tolerance * scale) { return ::testing::AssertionSuccess(); }
		return ::testing::AssertionFailure() << std::format(
			"The relative difference between {} and {} is {:.3e}, which exceeds {} = {:.3e}, where\n{} = {:.17g} and\n{} = {:.17g}",
			a_expression, b_expression, std::abs(a - b) / scale, tolerance_expression, tolerance, a_expression, a, b_expression, b
		);
	}
}

/// Expects the relative difference of the doubles \p a and \p b to be at most \p tolerance (see `dm::testing::relative_near`)
#define EXPECT_REL_NEAR(a, b, tolerance) EXPECT_PRED_FORMAT3(::dm::testing::relative_near, a, b, tolerance)
/// Asserts that the relative difference of the doubles \p a and \p b is at most \p tolerance (see `dm::testing::relative_near`)
#define ASSERT_REL_NEAR(a, b, tolerance) ASSERT_PRED_FORMAT3(::dm::testing::relative_near, a, b, tolerance)

#endif
