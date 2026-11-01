#ifndef DM_COMMON_CROSS_SECTION_FUNCTIONS_H
#define DM_COMMON_CROSS_SECTION_FUNCTIONS_H

#include <numbers>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/trf_kinematics.hpp"

namespace dm {

namespace cross_section_functions {
	/// Computes the cross-section prefactor `100 pi M_W^4 / (2 pi (Q^2 + M_W^2)^2 M E)` when `use_modified_prefactor == true` and `G_F^2 M_W^4 / (2 pi (Q^2 + M_W^2)^2)` otherwise
	constexpr double prefactor(const trf_kinematics &kinematics, const bool use_modified_prefactor) noexcept {
		if (use_modified_prefactor) {
			constexpr double numerator = 50.0 * math::power<4>(constants::particles::W.mass);
			const double denominator = math::power<2>(kinematics.Q2 + math::power<2>(constants::particles::W.mass)) * kinematics.target_mass * kinematics.E_beam;
			return numerator / denominator;
		} else {
			constexpr double numerator = math::power<2>(constants::fermi_coupling) * math::power<4>(constants::particles::W.mass);
			const double denominator = 2.0 * std::numbers::pi * math::power<2>(kinematics.Q2 + math::power<2>(constants::particles::W.mass));

			return numerator / denominator;
		}
	}
	/// Computes the Jacobian of the transformation `(x, Q^2) -> (x, y)`
	constexpr double xy_jacobian(const trf_kinematics &kinematics, const scattering_process &process) noexcept {
		return (kinematics.s - math::power<2>(process.target.mass) - math::power<2>(process.projectile.mass)) * kinematics.x;
	}
}

}

#endif