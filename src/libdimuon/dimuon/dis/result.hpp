#ifndef DM_DIS_RESULT_H
#define DM_DIS_RESULT_H

#include <flat_set>
#include <vector>

#include "dimuon/common/perturbative_quantity.hpp"
#include "dimuon/dis/physical_channel.hpp"

namespace dm {

namespace dis {
	/// A DIS cross section (or structure function) and its decomposition into the contributions of each physical channel (see `dis::evaluator::evaluate_observable`)
	struct result {
		/// Total value, the sum of the channel decomposition
		perturbative_quantity value;
		/// The physical channels of the decomposition (see `dis::physical_channel`), in the order of `channel_decomposition`
		std::flat_set<dis::physical_channel> physical_channels;
		/// The contribution of each physical channel, which holds the values of all perturbative orders
		std::vector<perturbative_quantity> channel_decomposition;

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
