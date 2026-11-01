#ifndef DM_DECAY_INSTANCE_H
#define DM_DECAY_INSTANCE_H

#include <algorithm>
#include <limits>
#include <type_traits>

#include "dimuon/decay/decay_function.hpp"

#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/common/particle.hpp"

namespace dm {

/// Encapsulation of the decay of a resonance, with an associated decay function of type \p decay_function_t
template <is_decay_function decay_function_t>
struct decay_instance {
	/// Resonance particle (which initiates the decay)
	particle resonance;
	/// Target hadron of the entire scattering process. 
	/// This is only needed for computing kinematical variables; the decay itself doesn't depend on the target.
	const particle hadron;
	/// Callable implementing the decay function
	const decay_function_t decay_function;

	/// Minimum lepton momentum cut-off (GeV), typically 5
	double lepton_momentum_min;
	/// Minimum `z_min` cut-off (typically 0; see `z_min()`)
	double z_min_cutoff;

	constexpr explicit decay_instance(
		particle resonance_, 
		particle hadron_, 
		decay_function_t decay_function_, 
		double lepton_momentum_min_, 
		double z_min_cutoff_ = 0.0) noexcept(std::is_nothrow_move_constructible_v<decay_function_t>)
	: resonance{std::move(resonance_)}, 
	hadron{std::move(hadron_)}, 
	decay_function{std::move(decay_function_)}, 
	lepton_momentum_min{lepton_momentum_min_}, 
	z_min_cutoff{z_min_cutoff_} {}

	constexpr double operator()(const double x, const double z, const double Q2) const noexcept(noexcept(decay_function(x, z, Q2, hadron))) {
		if (x != prev_x || z != prev_z || Q2 != prev_Q2) {
			prev_x = x;
			prev_z = z;
			prev_Q2 = Q2;

			prev_value = decay_function(x, z, Q2, hadron);
		}
		return prev_value;
	}

	/// Calculates `z_min` (i.e. the lower bound of the `z` integral)
	constexpr double z_min(const trf_kinematics &kinematics) const noexcept {
		return std::max({
			lepton_momentum_min / (kinematics.y * kinematics.E_beam), 
			resonance.mass / (kinematics.y * kinematics.E_beam),
			z_min_cutoff
		});
	}

private:
	mutable double prev_x = std::numeric_limits<double>::quiet_NaN();
	mutable double prev_z = std::numeric_limits<double>::quiet_NaN();
	mutable double prev_Q2 = std::numeric_limits<double>::quiet_NaN();
	mutable double prev_value = 0.0;
};

}

#endif