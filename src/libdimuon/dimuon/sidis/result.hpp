#ifndef DM_SIDIS_RESULT_H
#define DM_SIDIS_RESULT_H

#include <flat_set>

#include "dimuon/common/perturbative_quantity.hpp"
#include "dimuon/sidis/physical_channel.hpp"
#include "dimuon/utility/tensor.hpp"

namespace dm {

namespace sidis {
	/// A SIDIS cross section and its decomposition into the contributions of each pair of a row of the PDF values and a row of the FF values
	/// and each physical channel (see `sidis::evaluator::evaluate_observable`). The rows are the PDF and FF interpolation kernels in a kernel grid
	/// (with `vector_kernel_interface` in place of the PDF and the FFs), and otherwise a single row each, the PDF and the FFs.
	struct result {
		/// Total value
		perturbative_quantity value;
		/// The physical channels of the decomposition (see `sidis::physical_channel`), in the order of the last index of `channel_decomposition`
		std::flat_set<sidis::physical_channel> physical_channels;
		/// The contribution of each (PDF row, FF row, physical channel), which holds the values of all perturbative orders
		tensor<perturbative_quantity, 3> channel_decomposition;

		result &operator*=(const double rhs) noexcept {
			value *= rhs;
			for (perturbative_quantity &channel_value : channel_decomposition) {
				channel_value *= rhs;
			}
			return *this;
		}

		/// Product of the result and \p rhs (see `operator*=`)
		friend result operator*(result lhs, const double rhs) noexcept { lhs *= rhs; return lhs; }
		/// Product of \p lhs and the result (see `operator*=`)
		friend result operator*(const double lhs, result rhs) noexcept { rhs *= lhs; return rhs; }
	};
}

}

#endif
