#ifndef DM_COMMON_PARTICLE_H
#define DM_COMMON_PARTICLE_H

#include <iostream>
#include <limits>

namespace dm {

/// Container for the following particle data: mass, lifetime
struct particle {
	/// Mass of the particle (GeV)
	double mass = 0.0;
	/// Lifetime of the particle (ps). Infinite for a stable particle.
	double lifetime = std::numeric_limits<double>::infinity();
	/// Returns whether the particle is stable or not
	bool is_stable() const noexcept { return lifetime == std::numeric_limits<double>::infinity(); }

	bool operator==(const particle &) const = default;

	friend std::ostream& operator<<(std::ostream &os, const particle &o) {
		if (o.is_stable()) {
			return os << "m = " << o.mass;
		} else {
			return os << "m = " << o.mass << " | " << "tau = " << o.lifetime;
		}
	}
};

}

#endif