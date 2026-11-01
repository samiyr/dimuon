#ifndef DM_SIDIS_COEFFICIENT_CHANNEL_H
#define DM_SIDIS_COEFFICIENT_CHANNEL_H

#include <algorithm>
#include <compare> // IWYU pragma: keep (needed by the defaulted operator<=>)
#include <utility>

#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

namespace sidis {
	/// The properties of a channel that determine the coefficient function values.
	/// See the section "Channels" in `numerical_evaluation.md` for an explanation.
	struct coefficient_channel {
		/// Perturbative order of the channel
		perturbative_order pto;
		/// Type of the channel up to charge conjugation: the smaller of the type and the type of the charge-conjugate channel (see `partonic_channel::conjugate_type`).
		/// Charge-conjugate channels have the same coefficient functions (the sign of the antiquark terms in F3 is applied with the coupling).
		partonic_channel::channel_type c_symmetric_type;
		/// Mass scheme of the channel
		mass_scheme_t mass_scheme;
		/// Mass of the fragmenting parton
		double fragmentation_mass;
		/// Total (linear) mass of the final state (see `partonic_channel::final_state_mass`)
		double final_state_mass;

		/// Constructs the coefficient channel of \p channel
		explicit coefficient_channel(const partonic_channel &channel) noexcept :
			pto{channel.pto},
			c_symmetric_type{std::min(channel.type, partonic_channel::conjugate_type(channel.type))},
			mass_scheme{channel.mass_scheme},
			fragmentation_mass{channel.fragmentation_line.mass},
			final_state_mass{channel.final_state_mass()} {}

		/// Default the comparison operators to get member-wise ordering
		auto operator<=>(const coefficient_channel &) const = default;

		/// Decide whether subtraction terms (in a GM-VFNS) should be included based on the given (squared) energy scale.
		/// Subtraction terms are included if both of the following hold:
		/// - the mass scheme is a subtraction scheme
		/// - \p is above the (square of the) mass of the fragmenting parton.
		bool include_subtraction_term(const double scale) const noexcept {
			return math::power<2>(fragmentation_mass) < scale 
				&& (mass_scheme == mass_scheme_t::sacot_subtracted
					|| mass_scheme == mass_scheme_t::sacot_subtracted_modified
				);
		}

		/// Determine whether the mass scheme of this channel is a "modified" scheme.
		/// In a modified scheme, the F2 structure function is proportional to the Bjorken-x,
		/// as opposed to the slow-rescaling variable.
		bool is_modified_mass_scheme() const noexcept {
			return mass_scheme == mass_scheme_t::slow_rescaling_modified
				|| mass_scheme == mass_scheme_t::sacot_unsubtracted_modified
				|| mass_scheme == mass_scheme_t::sacot_subtracted_modified;
		}

		/// Computes the factor `1 + m^2 / Q^2`, where `m` is the final-state mass `final_state_mass`
		double chi_factor(const double Q2) const noexcept {
			return (1.0 + math::power<2>(final_state_mass) / Q2);
		}

		/// Gives the kinematical momentum fraction of this channel by taking potential mass effects into account
		double chi(const double x, const double Q2) const noexcept {
			return x * chi_factor(Q2);
		}
		/// Returns the appropriate PDF momentum fraction variable.
		/// In the ZM scheme, gives just `x`. In the SR and SACOT schemes, gives the slow-rescaling variable (i.e. `chi()`).
		double pdf_momentum_fraction(const double x, const double Q2) const noexcept {
			switch(mass_scheme) {
			using enum mass_scheme_t;
			case zero_mass:
				return x;

			case slow_rescaling:
			case sacot_unsubtracted:
			case sacot_subtracted:
			case slow_rescaling_modified:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return chi(x, Q2);
			}
			std::unreachable();
		}

		/// Returns true if both the initial and fragmenting partons are quarks
		bool is_quark_to_quark() const noexcept {
			return !is_quark_to_gluon() && !is_gluon_to_quark() && !is_gluon_to_gluon();
		}
		/// Returns true if the initial parton is a quark and the fragmenting parton is a gluon
		bool is_quark_to_gluon() const noexcept {
			using enum partonic_channel::channel_type;
			return c_symmetric_type == g_alpha || c_symmetric_type == g_betabar;
		}
		/// Returns true if the initial parton is a gluon and the fragmenting parton is a quark
		bool is_gluon_to_quark() const noexcept {
			using enum partonic_channel::channel_type;
			return c_symmetric_type == beta_g || c_symmetric_type == alphabar_g;
		}
		/// Returns true if both the initial and fragmenting partons are gluons
		bool is_gluon_to_gluon() const noexcept {
			return c_symmetric_type == partonic_channel::channel_type::gg;
		}
	};
}

}

#endif
