#ifndef DM_COUPLING_ELECTROMAGNETIC_COUPLING_H
#define DM_COUPLING_ELECTROMAGNETIC_COUPLING_H

#include "dimuon/common/flavor.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

struct electromagnetic_coupling {
	double operator()(const flavor::type flavor, const flavor::type other) const noexcept {
		if (flavor::abs(flavor) != flavor::abs(other)) { return 0.0; }
		const double charge = [flavor]() {
			if (flavor::is_upper_flavor(flavor)) { return 2.0 / 3.0; }
			if (flavor::is_lower_flavor(flavor)) { return -1.0 / 3.0; }
			return 0.0;
		}();

		return 0.5 * math::power<2>(charge);
	}
};

}

#endif