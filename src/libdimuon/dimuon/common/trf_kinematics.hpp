#ifndef DM_COMMON_TRF_KINEMATICS_H
#define DM_COMMON_TRF_KINEMATICS_H

#include "dimuon/utility/math.hpp"

namespace dm {

/// A structure encapsulating target rest-frame kinematics
struct trf_kinematics {
	/// Bjorken-x (momentum fraction), in the range [0, 1]
	const double x;
	/// Inelasticity, in the range [0, 1]
	const double y;
	/// Beam energy (GeV)
	const double E_beam;

	/// Target (i.e. proton or nucleus) mass (GeV)
	const double target_mass;
	/// Projectile (i.e. neutrino) mass (GeV)
	const double projectile_mass;

	/// Mandelstam s (GeV^2)
	const double s;
	/// Virtuality (GeV^2)
	const double Q2;

	/// @brief Constructs the kinematics from x, y and the beam energy
	/// @param[in] x Bjorken-x, must be in the range [0, 1]
	/// @param[in] y Inelasticity, must be in the range [0, 1]
	/// @param[in] E_beam Beam energy (GeV)
	/// @param[in] target_mass Target (i.e. proton or nucleus) mass (GeV)
	/// @param[in] projectile_mass Projectile (i.e. neutrino) mass (GeV)
	constexpr static trf_kinematics y_E_beam(
		const double x, const double y, const double E_beam, 
		const double target_mass, const double projectile_mass) noexcept {

		const double s = trf_kinematics::s_from_beam_energy(E_beam, target_mass, projectile_mass);
		const double Q2 = trf_kinematics::Q2_from_y(x, y, E_beam, target_mass);

		return trf_kinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	/// @brief Constructs the kinematics from x, Q^2 and s
	/// @param[in] x Bjorken-x, must be in the range [0, 1]
	/// @param[in] Q2 Virtuality (GeV^2)
	/// @param[in] s Mandelstam s (GeV^2)
	/// @param[in] target_mass Target (i.e. proton or nucleus) mass (GeV)
	/// @param[in] projectile_mass Projectile (i.e. neutrino) mass (GeV)
	constexpr static trf_kinematics Q2_s(
		const double x, const double Q2, const double s, 
		const double target_mass, const double projectile_mass) noexcept {

		const double y = y_from_Q2(x, Q2, s, target_mass, projectile_mass);
		const double E_beam = beam_energy_from_s(s, target_mass, projectile_mass);

		return trf_kinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	/// @brief Constructs the kinematics from x, Q^2 and sqrt(s)
	/// @param[in] x Bjorken-x, must be in the range [0, 1]
	/// @param[in] Q2 Virtuality (GeV^2)
	/// @param[in] sqrt_s Square root of Mandelstam s (GeV)
	/// @param[in] target_mass Target (i.e. proton or nucleus) mass (GeV)
	/// @param[in] projectile_mass Projectile (i.e. neutrino) mass (GeV)
	constexpr static trf_kinematics Q2_sqrt_s(
		const double x, const double Q2, const double sqrt_s, 
		const double target_mass, const double projectile_mass) noexcept {
		return trf_kinematics::Q2_s(x, Q2, sqrt_s * sqrt_s, target_mass, projectile_mass);
	}

	/// Checks whether the kinematical variables are in their physically valid ranges
	constexpr bool is_valid() const noexcept {
		return (y > 0.0 && y <= 1.0 && x > 0.0 && x <= 1.0 && E_beam > 0.0 && target_mass >= 0.0 && projectile_mass >= 0.0 && s >= 0.0 && Q2 >= 0.0);
	}

private:
	/// Conversion from beam energy to Mandelstam s
	constexpr static double s_from_beam_energy(const double E_beam, const double target_mass, const double projectile_mass) noexcept {
		return math::power<2>(target_mass) + math::power<2>(projectile_mass) + 2 * target_mass * E_beam;
	}
	/// Conversion from Mandelstam s to beam energy
	constexpr static double beam_energy_from_s(const double s, const double target_mass, const double projectile_mass) noexcept {
		return (s - math::power<2>(target_mass) - math::power<2>(projectile_mass)) / (2 * target_mass);
	}
	/// Conversion from y to Q^2
	constexpr static double Q2_from_y(const double x, const double y, const double E_beam, const double target_mass) noexcept {
		return 2 * target_mass * E_beam * x * y;
	}
	/// Conversion from Q^2 to y
	constexpr static double y_from_Q2(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass) noexcept {
		return Q2 / ((s - math::power<2>(target_mass) - math::power<2>(projectile_mass)) * x);
	}
};

}

#endif