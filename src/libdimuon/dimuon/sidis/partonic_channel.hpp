#ifndef DM_SIDIS_PARTONIC_CHANNEL_H
#define DM_SIDIS_PARTONIC_CHANNEL_H

#include <iostream>
#include <string>
#include <string_view>
#include <format>
#include <numeric>
#include <utility>
#include <ranges>
#include <algorithm>

#include "dimuon/common/flavor.hpp"
#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/common/parton.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/utility/collections.hpp"

namespace dm {

namespace sidis {
	/// Representation of an individual partonic SIDIS channel
	struct partonic_channel {
		/// Types of SIDIS channels, based on the notation in 2506.19926 [hep-ph]:
		/// alpha: lower-type quark (d, s, b)
		/// beta: upper-type quark (u, c)
		/// q: quark of arbitrary flavor
		/// bar: antiquark
		enum class channel_type {
			/// alpha --> beta
			beta_alpha,
			/// q --> q (flavor-conserving)
			qq_Fcon,
			/// alpha --> alpha (flavor-conserving)
			alpha_alpha_Fcon,
			/// alpha --> anti-beta
			betabar_alpha,
			/// alpha --> anti-alpha (flavor-conserving)
			alphabar_alpha_Fcon,
			/// alpha --> q
			qprime_alpha,
			/// q --> beta
			betaprime_q,
			/// alpha --> anti-q
			qbarprime_alpha,
			/// q --> anti-alpha
			alphabarprime_q,
			/// alpha --> gluon
			g_alpha,

			/// anti-beta --> anti-alpha [conjugate to `beta_alpha`]
			alphabar_betabar,
			/// anti-q --> anti-q (flavor-conserving) [conjugate to `qq_Fcon`]
			qbar_qbar_Fcon,
			/// anti-beta --> anti-beta (flavor-conserving) [conjugate to `alpha_alpha_Fcon`]
			betabar_betabar_Fcon,
			/// anti-beta --> alpha [conjugate to `betabar_alpha`]
			alpha_betabar,
			/// anti-beta --> beta (flavor-conserving) [conjugate to `alphabar_alpha_Fcon`]
			beta_betabar_Fcon,
			/// anti-beta --> q [conjugate to `qprime_alpha`]
			qprime_betabar,
			/// anti-q --> beta [conjugate to `betaprime_q`]
			betaprime_qbar,
			/// anti-beta --> anti-q [conjugate to `qbarprime_alpha`]
			qbarprime_betabar,
			/// anti-q --> anti-alpha [conjugate to `alphabarprime_q`]
			alphabarprime_qbar,
			/// anti-beta --> gluon [conjugate to `g_alpha`]
			g_betabar,

			/// gluon --> beta
			beta_g,
			/// gluon --> anti-alpha [conjugate to `beta_g`]
			alphabar_g,

			/// gluon --> gluon
			gg
		};
		perturbative_order pto;

		/// Type of channel
		channel_type type;
		/// The initial-state parton connected to the PDF
		parton initial_line;
		/// The final-state parton connected to the FF
		parton fragmentation_line;
		/// Array of final-state partons that do not fragment.
		/// Used in the case where there is a quark which
		/// is not connected to the FF. At NLO, this means
		/// the diagrams where a gluon is fragmenting or
		/// the diagrams with an initial-state gluon and
		/// two outgoing quarks. 
		std::vector<parton> nonfragmentation_lines;
		/// Mass scheme of the channel
		mass_scheme_t mass_scheme;
		/// Partons that "really" fragment, used when fragmenting a
		/// non-charm parton to get the correct number of final-state charm quarks
		std::vector<parton> true_fragmentation_lines;

		/// @brief Default constructor which initializes each member variable
		/// @param[in] pto_ Perturbative order of the channel
		/// @param[in] type_ Type of channel
		/// @param[in] initial_line_ The initial-state parton connected to the PDF
		/// @param[in] fragmentation_line_ The final-state parton connected to the FF
		/// @param[in] nonfragmentation_lines_ Array of final-state partons that do not fragment
		/// @param[in] mass_scheme_ Mass scheme of the channel
		/// @param[in] true_fragmentation_lines_ Partons that "really" fragment (see member field description)
		explicit partonic_channel(
			perturbative_order pto_,
			channel_type type_,
			parton initial_line_, 
			parton fragmentation_line_, 
			std::vector<parton> nonfragmentation_lines_, 
			mass_scheme_t mass_scheme_,
			std::vector<parton> true_fragmentation_lines_
		) noexcept :
		pto{pto_},
		type{type_},
		initial_line{std::move(initial_line_)}, 
		fragmentation_line{std::move(fragmentation_line_)}, 
		nonfragmentation_lines{std::move(nonfragmentation_lines_)},
		mass_scheme{mass_scheme_},
		true_fragmentation_lines{std::move(true_fragmentation_lines_)} {}

		/// Default the equality operator to get the compiler-generated member-wise equality comparison
		bool operator==(const partonic_channel &) const = default;

		/// Returns the type of the charge-conjugate channel of a channel of type \p type
		static channel_type conjugate_type(const channel_type type) noexcept {
			switch (type) {
				using enum channel_type;
				case beta_alpha: 			return alphabar_betabar;
				case qq_Fcon: 				return qbar_qbar_Fcon;
				case alpha_alpha_Fcon: 		return betabar_betabar_Fcon;
				case betabar_alpha: 		return alpha_betabar;
				case alphabar_alpha_Fcon: 	return beta_betabar_Fcon;
				case qprime_alpha: 			return qprime_betabar;
				case betaprime_q: 			return betaprime_qbar;
				case qbarprime_alpha: 		return qbarprime_betabar;
				case alphabarprime_q: 		return alphabarprime_qbar;
				case g_alpha: 				return g_betabar;

				case alphabar_betabar: 		return beta_alpha;
				case qbar_qbar_Fcon: 		return qq_Fcon;
				case betabar_betabar_Fcon: 	return alpha_alpha_Fcon;
				case alpha_betabar: 		return betabar_alpha;
				case beta_betabar_Fcon: 	return alphabar_alpha_Fcon;
				case qprime_betabar: 		return qprime_alpha;
				case betaprime_qbar: 		return betaprime_q;
				case qbarprime_betabar: 	return qbarprime_alpha;
				case alphabarprime_qbar: 	return alphabarprime_q;
				case g_betabar: 			return g_alpha;

				case beta_g: 				return alphabar_g;
				case alphabar_g: 			return beta_g;

				case gg: 					return gg;
			}
			std::unreachable();
		}

		/// Computes the total (linear) mass of the final state
		double final_state_mass() const noexcept {
			const auto accumulator = [](double acc, const parton p) { return acc + p.mass; };

			const double fragmentation_mass = std::accumulate(true_fragmentation_lines.begin(), true_fragmentation_lines.end(), 0.0, accumulator);
			const double nonfragmentation_mass = std::accumulate(nonfragmentation_lines.begin(), nonfragmentation_lines.end(), 0.0, accumulator);

			return fragmentation_mass + nonfragmentation_mass;
		}

		/// Check that the electric charge is conserved in the channel. The sign of the W boson is determined from the \p process argument.
		bool check_charge_conservation(const scattering_process &process) const noexcept {
			const int W = 3 * process.W_sign();
			const int initial_charge = W + initial_line.charge();
			const int final_charge = fragmentation_line.charge() + std::ranges::fold_left(nonfragmentation_lines | std::views::transform([](const parton &p) { return p.charge(); }), 0, std::plus{});
			return initial_charge == final_charge;
		}

		/// Returns true if there is further fragmentation beyond the `fragmentation_line` (see `true_fragmentation_lines`)
		bool further_fragmentation() const noexcept {
			return !(true_fragmentation_lines.size() == 1 && true_fragmentation_lines[0] == fragmentation_line);
		}
		/// Returns a pair of partons that are connected to the gauge boson
		std::pair<parton, parton> get_coupling_quarks() const noexcept {
			switch(type) {
				case channel_type::beta_alpha: 			
				case channel_type::alphabar_betabar:
					return std::make_pair(initial_line, fragmentation_line);

				case channel_type::qq_Fcon:
				case channel_type::qbar_qbar_Fcon:				
					return std::make_pair(nonfragmentation_lines[0], nonfragmentation_lines[1]);

				case channel_type::alpha_alpha_Fcon:
				case channel_type::betabar_betabar_Fcon:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::betabar_alpha:
				case channel_type::alpha_betabar:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::alphabar_alpha_Fcon:
				case channel_type::beta_betabar_Fcon:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::qprime_alpha:
				case channel_type::qprime_betabar:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::betaprime_q:
				case channel_type::betaprime_qbar:
					return std::make_pair(fragmentation_line, nonfragmentation_lines[0]);

				case channel_type::qbarprime_alpha:
				case channel_type::qbarprime_betabar:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::alphabarprime_q:
				case channel_type::alphabarprime_qbar:
					return std::make_pair(fragmentation_line, nonfragmentation_lines[0]);

				case channel_type::g_alpha:
				case channel_type::g_betabar:
					return std::make_pair(initial_line, nonfragmentation_lines[0]);

				case channel_type::beta_g:
				case channel_type::alphabar_g:
					return std::make_pair(fragmentation_line, nonfragmentation_lines[0]);

				case channel_type::gg:
					return std::make_pair(nonfragmentation_lines[0], nonfragmentation_lines[1]);
			}
			std::unreachable();
		}

		/// Functor implementing hashing for `partonic_channel`
		struct hasher {
			std::size_t operator()(const partonic_channel &channel) const noexcept {
				using namespace collections;
				
				const auto hash_parton = [](std::size_t &seed, const parton &p) {
					hash_combine(seed, std::to_underlying(p.flavor));
					hash_combine(seed, p.mass);
				};

				std::size_t seed = 0;

				hash_combine(seed, std::to_underlying(channel.pto));
				hash_combine(seed, std::to_underlying(channel.type));

				hash_parton(seed, channel.initial_line);
				hash_parton(seed, channel.fragmentation_line);

				for (const parton &p : channel.nonfragmentation_lines) {
					hash_parton(seed, p);
				}

				hash_combine(seed, std::to_underlying(channel.mass_scheme));

				for (const parton &p : channel.true_fragmentation_lines) {
					hash_parton(seed, p);
				}

				return seed;
			}
		};
	};
}

}

/// Formats a SIDIS channel as its type, its mass scheme in parentheses, and its partons with their masses, separated by tabs.
/// The fragmenting parton is marked with `*`, followed by the partons that really fragment if they differ from it (see `further_fragmentation`).
template <> struct std::formatter<dm::sidis::partonic_channel> : std::formatter<std::string_view> {
	auto format(const dm::sidis::partonic_channel &channel, auto &context) const {
		const std::string_view channel_name = [&]() {
			switch (channel.type) {
				using enum dm::sidis::partonic_channel::channel_type;
				case beta_alpha: 				return "    α --> β    ";
				case qq_Fcon: 					return "    q --> q    ";
				case alpha_alpha_Fcon: 			return "    α --> α    ";
				case betabar_alpha: 			return "    α --> βbar ";
				case alphabar_alpha_Fcon: 		return "    α --> αbar ";
				case qprime_alpha: 				return "   α  --> q'   ";
				case betaprime_q: 				return "   q  --> β'   ";
				case qbarprime_alpha: 			return "    α --> qbar'";
				case alphabarprime_q: 			return "    q --> αbar'";
				case g_alpha: 					return "    α --> g    ";

				case alphabar_betabar: 			return " βbar --> αbar ";
				case qbar_qbar_Fcon: 			return " qbar --> qbar ";
				case betabar_betabar_Fcon: 		return " βbar --> βbar ";
				case alpha_betabar: 			return " βbar --> α    ";
				case beta_betabar_Fcon: 		return " βbar --> β    ";
				case qprime_betabar: 			return " βbar --> q'   ";
				case betaprime_qbar: 			return " qbar --> β    ";
				case qbarprime_betabar: 		return " βbar --> qbar'";
				case alphabarprime_qbar: 		return " qbar --> αbar'";
				case g_betabar: 				return " βbar --> g    ";

				case beta_g: 					return "    g --> β    ";
				case alphabar_g: 				return "    g --> αbar ";
				case gg: 						return "    g --> g    ";
			}
			std::unreachable();
		}();
		std::string channel_flow = std::format("{} [m = {}] --> ", channel.initial_line.flavor, channel.initial_line.mass);

		for (const dm::parton &p : channel.nonfragmentation_lines) {
			channel_flow += std::format("{} [m = {}] + ", p.flavor, p.mass);
		}

		channel_flow += std::format("{} [m = {}]*", channel.fragmentation_line.flavor, channel.fragmentation_line.mass);

		if (channel.further_fragmentation()) {
			channel_flow += " { -->";
			for (const dm::parton &p : channel.true_fragmentation_lines) {
				channel_flow += std::format(" {} [m = {}]*", p.flavor, p.mass);
			}
			channel_flow += " }";
		}

		return std::formatter<std::string_view>::format(std::format("{}\t({})\t{}", channel_name, channel.mass_scheme, channel_flow), context);
	}
};

namespace dm {

namespace sidis {
	/// Prints the channel as it is formatted (see `std::formatter<dm::sidis::partonic_channel>`)
	inline std::ostream &operator<<(std::ostream &os, const partonic_channel &channel) {
		return os << std::format("{}", channel);
	}
}

}

#endif