#ifndef DM_COMMON_SPLITTING_H
#define DM_COMMON_SPLITTING_H

#include "dimuon/utility/math.hpp"

namespace dm {

namespace splitting {
	constexpr double Pqg(const double x) noexcept {
		return 0.5 * (math::power<2>(x) + math::power<2>(1.0 - x));
	}
}

}

#endif