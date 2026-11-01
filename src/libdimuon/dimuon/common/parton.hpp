#ifndef DM_COMMON_PARTON_H
#define DM_COMMON_PARTON_H

#include "dimuon/common/flavor.hpp"

namespace dm {

/// Representation of an individual parton
struct parton {
	/// Flavor of the parton, given by its PID (see `flavor` for more)
	flavor::type flavor;
	/// Mass of the parton (GeV)
	double mass = 0.0;
	/// Returns the charge of the parton, encoded by `flavor::charge`
	int charge() const noexcept {
		return flavor::charge(flavor);
	}

	bool operator==(const parton &) const = default;

	/// Returns the conjugate of the parton \p p (i.e. a parton of the same mass but conjugate flavor).
	/// The conjugate flavor converts between quarks and antiquarks; see `flavor::conjugate_flavor`.
	static parton conjugate_parton(const parton &p) noexcept {
		return parton { flavor::conjugate_flavor(p.flavor), p.mass };
	}
};

}

#endif