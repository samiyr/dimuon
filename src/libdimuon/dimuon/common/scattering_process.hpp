#ifndef DM_COMMON_SCATTERING_PROCESS_H
#define DM_COMMON_SCATTERING_PROCESS_H

#include <iostream>
#include <utility>

#include "dimuon/common/particle.hpp"

namespace dm {

/// A structure encapsulating a reaction process
struct scattering_process {
	/// Type of the process, either neutrino or antineutrino scattering
	enum class process_type {
		/// neutrino + target --> lepton + X
		neutrino,
		/// antineutrino + target --> antilepton + X
		antineutrino,
		/// lepton + target --> lepton + X
		electromagnetic
	};

	/// Process type, either neutrino or antineutrino scattering
	process_type type;
	/// Target particle (i.e. the proton or a nucleus)
	particle target;
	/// Projectile particle (i.e. neutrino or antineutrino)
	particle projectile;

	/// Gives the sign of the process (+1 for neutrino scattering, -1 for antineutrino scattering, 0 for electromagnetic scattering) 
	/// used in the cross section formula (see `cross_section_observable`)
	constexpr int W_sign() const noexcept {
		switch (type) {
			case process_type::neutrino: return 1;
			case process_type::antineutrino: return -1;
			case process_type::electromagnetic: return 0;
		}

		std::unreachable();
	}

	/// Checks whether the sign of the process is positive or not; see scattering_process::W_sign
	constexpr bool positive_W() const noexcept {
		return W_sign() == 1;
	}

	friend std::ostream& operator<<(std::ostream &os, const scattering_process &o) {
		std::string type_string;
		switch (o.type) {
		case process_type::neutrino:
			type_string = "neutrino to lepton";
			break;
		case process_type::antineutrino:
			type_string = "antineutrino to antilepton";
			break;
		case process_type::electromagnetic:
			type_string = "electromagnetic";
			break;
		}

		return os << type_string << " [target = (" << o.target << ") | projectile = (" << o.projectile << ")]";
	}
};

}

#endif