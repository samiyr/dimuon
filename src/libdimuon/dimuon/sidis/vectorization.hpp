#ifndef DM_SIDIS_VECTORIZATION_H
#define DM_SIDIS_VECTORIZATION_H

#include <cstddef>
#include <utility>

#include "dimuon/utility/vectorization.hpp"

namespace dm {

/// Helper functions for evaluating the SIDIS coefficient functions for all products of a PDF value and an FF value at once, with one product in each lane of a vector
namespace sidis::vectorization {
	/// Number of distinct PDF values needed in the coefficient functions: the values at `momentum_fraction` and `momentum_fraction / xi` (`pdf` and `pdf_hat`)
	inline constexpr std::size_t pdf_value_count = 2;
	/// Number of distinct FF values needed in the coefficient functions: the values at `z`, `z / zeta`, and `z / zeta_p` (`ff`, `ff_hat`, and `ff_hat_p`)
	inline constexpr std::size_t ff_value_count = 3;
	/// Number of lanes of `distribution_product`, one for each product of a PDF value and an FF value
	inline constexpr std::size_t lane_count = pdf_value_count * ff_value_count;

	/// Type for holding the products of all necessary PDF and FF values at a single integration point. 
	/// See section "Vectorization and masking" in the `numerical_evaluation.md` document for more.
	using distribution_product = dm::vectorization::vec<lane_count>;

	/// The PDF and FF values needed in the coefficient functions (see the fields of the same names in `sidis::evaluation_parameters`)
	enum class distribution_value { pdf, pdf_hat, ff, ff_hat, ff_hat_p };

	/// Returns the unit value (masking vector) of \p value: one in the lanes of the products that contain it, and zero elsewhere.
	/// For the PDF value `a` this is the Kronecker product `v_a * 1` and for the FF value `b` it is `1 * u_b`, where `v_a` and `u_b` are unit vectors and `1` are vectors of ones.
	/// The element-wise product of the unit values of `a` and `b` is then `v_a * u_b`, which is one only in the lane of their product.
	/// See section "Vectorization and masking" in the `numerical_evaluation.md` document for more.
	inline distribution_product make_masking_vector(const distribution_value value) noexcept {
		using dm::vectorization::kronecker_product, dm::vectorization::unit_vector, dm::vectorization::ones;
		switch (value) {
			using enum distribution_value;
			case pdf: 		return kronecker_product(unit_vector<pdf_value_count>(0), ones<ff_value_count>());
			case pdf_hat: 	return kronecker_product(unit_vector<pdf_value_count>(1), ones<ff_value_count>());
			case ff: 		return kronecker_product(ones<pdf_value_count>(), unit_vector<ff_value_count>(0));
			case ff_hat: 	return kronecker_product(ones<pdf_value_count>(), unit_vector<ff_value_count>(1));
			case ff_hat_p: 	return kronecker_product(ones<pdf_value_count>(), unit_vector<ff_value_count>(2));
		}
		std::unreachable();
	}
}

}

#endif
