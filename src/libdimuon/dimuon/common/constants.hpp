#ifndef DM_COMMON_CONSTANTS_H
#define DM_COMMON_CONSTANTS_H

#include <cstddef>

#include "dimuon/common/particle.hpp"

#include "dimuon/utility/math.hpp"

namespace dm {

/// Provides access to certain physical/mathematical constants
namespace constants {
	/// QCD color factor (fundamental representation)
	inline constexpr double C_F = 4.0 / 3.0;
	/// QCD color factor (adjoint representation)
	inline constexpr double C_A = 3.0;
	/// QCD normalization of tr(t^a t^b) = T_R delta^ab in the fundamental representation
	inline constexpr double T_R = 0.5;
	/// Fermi coupling constant
	inline constexpr double fermi_coupling = 1.1663787e-5;

	/// The squared CKM matrix elements

	inline constexpr double V_ud_sq = math::power<2>(0.97435);
	inline constexpr double V_us_sq = math::power<2>(0.22500);
	inline constexpr double V_ub_sq = math::power<2>(0.00369);

	inline constexpr double V_cd_sq = math::power<2>(0.22486);
	inline constexpr double V_cs_sq = math::power<2>(0.97349);
	inline constexpr double V_cb_sq = math::power<2>(0.04182);

	inline constexpr double V_td_sq = math::power<2>(0.00857);
	inline constexpr double V_ts_sq = math::power<2>(0.04110);
	inline constexpr double V_tb_sq = math::power<2>(0.999112);

	/// Particle properties
	namespace particles {
		inline constexpr particle proton{0.938272};
		inline constexpr particle D0{1.86484, 0.4101};
		inline constexpr particle Dp{1.86962, 1.04};
		inline constexpr particle Ds{1.96847, 0.5};
		inline constexpr particle LambdaC{2.28646, 0.2};
		inline constexpr particle W{80.377, 4.7961631e-9};
		inline constexpr particle neutrino{0.0};
		inline constexpr particle muon{0.105658372};
		inline constexpr particle massless_muon{0.0};
	}

	/// Value of the Riemann zeta function zeta(3)
	inline constexpr double zeta3 = 1.20205690315959428539973816151;

	/// Total number of flavors, including the gluon and antiquark flavors
	inline constexpr std::size_t total_flavors = 13;

	/// Number of colors
	inline constexpr double NC = 3.0;
} 

}

#endif