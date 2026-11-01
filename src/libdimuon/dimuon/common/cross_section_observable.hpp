#ifndef DM_COMMON_CROSS_SECTION_OBSERVABLE_H
#define DM_COMMON_CROSS_SECTION_OBSERVABLE_H

#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/common/scattering_process.hpp"

namespace dm {

/// Functor type that calculates the charged-current cross section `d^2 sigma/dx dQ^2` from given F1, F2, and F3 structure function values 
struct cross_section_observable {
	double operator()(const trf_kinematics &kinematics, const scattering_process &process, const double f1, const double f2, const double f3) const noexcept {
		const double x = kinematics.x;
		const double y = kinematics.y;
		const double s = kinematics.s;
		const double M2 = math::power<2>(process.target.mass);

		const double term1 = x * math::power<2>(y);
		const double term2 = 1.0 - y - (x * y * M2) / (s - M2);
		const double term3 = y * (1 - 0.5 * y);

		const double result = (term1 * f1 + term2 * f2) / x + double(process.W_sign()) * term3 * f3;

		return result;
	}
};

}

#endif