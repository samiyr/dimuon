#ifndef DM_DIS_PARTONIC_CHANNEL_H
#define DM_DIS_PARTONIC_CHANNEL_H

#include <algorithm>
#include <format>
#include <functional>
#include <iostream>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "dimuon/common/flavor.hpp"
#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/common/parton.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/scattering_process.hpp"

namespace dm {

namespace dis {
	/// Representation of an individual partonic DIS channel
	struct partonic_channel {
		/// Types of DIS channels, in the notation of the SIDIS channels (see `sidis::partonic_channel::channel_type`), named as in neutrino scattering:
		/// alpha: lower-type quark (d, s, b)
		/// beta: upper-type quark (u, c)
		/// bar: antiquark
		enum class channel_type {
			/// alpha --> beta, at NLO with the emission of a gluon
			beta_alpha,
			/// anti-beta --> anti-alpha [conjugate to `beta_alpha`]
			alphabar_betabar,
			/// gluon --> beta + anti-alpha
			beta_g,
			/// gluon --> anti-alpha + beta [conjugate to `beta_g`]
			alphabar_g
		};
		/// Perturbative order of the channel
		perturbative_order pto;
		/// Type of channel
		channel_type type;
		/// The initial-state parton connected to the PDF
		parton initial_line;
		/// The final-state quarks: the quark into which the gauge boson turns the initial quark, or the quark and the antiquark of the gluon
		std::vector<parton> final_lines;
		/// Mass scheme of the channel
		mass_scheme_t mass_scheme;

		/// @brief Constructor which initializes each member variable
		/// @param[in] pto_ Perturbative order of the channel
		/// @param[in] type_ Type of channel
		/// @param[in] initial_line_ The initial-state parton connected to the PDF
		/// @param[in] final_lines_ The final-state quarks (see `final_lines`)
		/// @param[in] mass_scheme_ Mass scheme of the channel
		explicit partonic_channel(
			perturbative_order pto_,
			channel_type type_,
			parton initial_line_,
			std::vector<parton> final_lines_,
			mass_scheme_t mass_scheme_
		) noexcept :
		pto{pto_},
		type{type_},
		initial_line{std::move(initial_line_)},
		final_lines{std::move(final_lines_)},
		mass_scheme{mass_scheme_} {}

		/// Default the equality operator to get the compiler-generated member-wise equality comparison
		bool operator==(const partonic_channel &) const = default;

		/// Returns the type of the charge-conjugate channel of a channel of type \p type
		static channel_type conjugate_type(const channel_type type) noexcept {
			switch (type) {
				using enum channel_type;
				case beta_alpha: 		return alphabar_betabar;
				case alphabar_betabar: 	return beta_alpha;
				case beta_g: 			return alphabar_g;
				case alphabar_g: 		return beta_g;
			}
			std::unreachable();
		}

		/// Computes the total (linear) mass of the final state, which determines the mass dependence of the coefficient functions.
		/// At most one of the final-state quarks is assumed to be massive.
		double final_state_mass() const noexcept {
			return std::ranges::fold_left(final_lines | std::views::transform([](const parton &p) { return p.mass; }), 0.0, std::plus{});
		}

		/// Computes the total (linear) mass of the initial and final states, which determines the slow-rescaling variable (see `dis::coefficient_channel::chi`)
		double total_mass() const noexcept {
			return initial_line.mass + final_state_mass();
		}

		/// Check that the electric charge is conserved in the channel. The sign of the W boson is determined from the \p process argument.
		bool check_charge_conservation(const scattering_process &process) const noexcept {
			const int W = 3 * process.W_sign();
			const int initial_charge = W + initial_line.charge();
			const int final_charge = std::ranges::fold_left(final_lines | std::views::transform([](const parton &p) { return p.charge(); }), 0, std::plus{});
			return initial_charge == final_charge;
		}

		/// Returns a pair of partons that are connected to the gauge boson
		std::pair<parton, parton> get_coupling_quarks() const noexcept {
			switch (type) {
				using enum channel_type;
				case beta_alpha:
				case alphabar_betabar:
					return std::make_pair(initial_line, final_lines[0]);

				case beta_g:
				case alphabar_g:
					return std::make_pair(final_lines[0], final_lines[1]);
			}
			std::unreachable();
		}
	};
}

}

/// Formats a DIS channel as its type, its mass scheme in parentheses, and its partons with their masses, separated by tabs
template <> struct std::formatter<dm::dis::partonic_channel> : std::formatter<std::string_view> {
	auto format(const dm::dis::partonic_channel &channel, auto &context) const {
		const std::string_view channel_name = [&]() {
			switch (channel.type) {
				using enum dm::dis::partonic_channel::channel_type;
				case beta_alpha: 		return "    α --> β    ";
				case alphabar_betabar: 	return " βbar --> αbar ";
				case beta_g: 			return "    g --> β    ";
				case alphabar_g: 		return "    g --> αbar ";
			}
			std::unreachable();
		}();

		std::string channel_flow = std::format("{} [m = {}] --> ", channel.initial_line.flavor, channel.initial_line.mass);
		for (const auto &[index, p] : std::views::enumerate(channel.final_lines)) {
			channel_flow += std::format("{}{} [m = {}]", index == 0 ? "" : " + ", p.flavor, p.mass);
		}

		return std::formatter<std::string_view>::format(std::format("{}\t({})\t{}", channel_name, channel.mass_scheme, channel_flow), context);
	}
};

namespace dm {

namespace dis {
	/// Prints the channel as it is formatted (see `std::formatter<dm::dis::partonic_channel>`)
	inline std::ostream &operator<<(std::ostream &os, const partonic_channel &channel) {
		return os << std::format("{}", channel);
	}
}

}

#endif
