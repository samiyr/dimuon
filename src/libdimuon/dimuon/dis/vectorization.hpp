#ifndef DM_DIS_VECTORIZATION_H
#define DM_DIS_VECTORIZATION_H

#include <cstddef>
#include <utility>

#include "dimuon/utility/vectorization.hpp"

namespace dm {

/// Helper functions for evaluating the DIS coefficient functions for all PDF values at once, with one PDF value in each lane of a vector
namespace dis::vectorization {
	/// Number of distinct PDF values needed in the coefficient functions: the values at `momentum_fraction` and `momentum_fraction / xi` (`pdf` and `pdf_hat`)
	inline constexpr std::size_t pdf_value_count = 2;

	/// Type for holding the coefficients of all necessary PDF values at a single integration point, one lane per PDF value.
	/// See section "Vectorization and masking" in the `numerical_evaluation.md` document for the SIDIS version, of which this is the one without the FFs.
	using distribution_product = dm::vectorization::vec<pdf_value_count>;

	/// The PDF values needed in the coefficient functions (see the fields of the same names in `dis::evaluation_parameters`)
	enum class distribution_value { pdf, pdf_hat };

	/// Returns the unit value (masking vector) of \p value: one in the lane of the PDF value and zero in the other
	inline distribution_product make_masking_vector(const distribution_value value) noexcept {
		switch (value) {
			using enum distribution_value;
			case pdf: 		return dm::vectorization::unit_vector<pdf_value_count>(0);
			case pdf_hat: 	return dm::vectorization::unit_vector<pdf_value_count>(1);
		}
		std::unreachable();
	}
}

}

#endif
