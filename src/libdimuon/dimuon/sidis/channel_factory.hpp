#ifndef DM_SIDIS_CHANNEL_FACTORY_H
#define DM_SIDIS_CHANNEL_FACTORY_H

#include <vector>
#include <optional>
#include <ranges>
#include <unordered_map>

#include "dimuon/common/channel_configuration.hpp"
#include "dimuon/common/flavor.hpp"
#include "dimuon/common/parton_configuration.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/sidis/partonic_channel.hpp"
#include "dimuon/utility/collections.hpp"

namespace dm {

namespace sidis {
	/// Factory class for constructing SIDIS channels, with separation for massless (~ no charm) and massive (~ charm) channels
	struct channel_factory {
		using enum sidis::partonic_channel::channel_type;

		/// The configurations of the channel types: a channel type that is missing, or whose configuration is `std::nullopt`, is not included
		using channel_configuration_list = std::unordered_map<sidis::partonic_channel::channel_type, std::optional<channel_configuration>>;

		/// Returns the types of the channels at the perturbative order \p pto
		static std::vector<sidis::partonic_channel::channel_type> channel_types(const perturbative_order pto) {
			switch (pto) {
				case perturbative_order::lo:
					return { beta_alpha, alphabar_betabar };
				case perturbative_order::nlo:
					return { beta_alpha, alphabar_betabar, g_alpha, g_betabar, beta_g, alphabar_g };
				case perturbative_order::nnlo:
					return {
						beta_alpha, qq_Fcon, alpha_alpha_Fcon, betabar_alpha, alphabar_alpha_Fcon, qprime_alpha, betaprime_q, qbarprime_alpha, alphabarprime_q,
						alphabar_betabar, qbar_qbar_Fcon, betabar_betabar_Fcon, alpha_betabar, beta_betabar_Fcon, qprime_betabar, betaprime_qbar, qbarprime_betabar, alphabarprime_qbar,
						g_alpha, g_betabar,
						beta_g, alphabar_g,
						gg
					};
			}
			std::unreachable();
		}

		/// @brief Builds the list of channels of all types: the quark-to-quark, quark-to-gluon, gluon-to-quark, and gluon-to-gluon channels, in this order
		/// @param[in] pto Perturbative order
		/// @param[in] process The scattering process, used to determine the sign of the gauge boson
		/// @param[in] flavors Parton configuration, which determines the allowed flavors and the masses
		/// @param[in] massless A map of channel configurations for each SIDIS channel type, used in massless cases (~ no charm fragmentation).
		/// If `std::nullopt` is given as the value, or the channel type is missing, that channel is not included.
		/// @param[in] massive A map of channel configurations for each SIDIS channel type, used in massive cases (~ massive charm fragmentation).
		/// If `std::nullopt` is given as the value, or the channel type is missing, that channel is not included.
		static std::vector<sidis::partonic_channel> create_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			return collections::join(
				create_quark_to_quark_channels(pto, process, flavors, massless, massive),
				create_quark_to_gluon_channels(pto, process, flavors, massless, massive),
				create_gluon_to_quark_channels(pto, process, flavors, massless, massive),
				create_gluon_to_gluon_channels(pto, process, flavors, massless, massive)
			);
		}

		/// Builds the list of all channels of the perturbative order \p pto (see `channel_types`), with the configuration \p massless
		/// for the massless channels and \p massive for the massive channels. If a configuration is `std::nullopt`, those channels are not included.
		static std::vector<sidis::partonic_channel> create_all_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const std::optional<channel_configuration> &massless,
			const std::optional<channel_configuration> &massive
		) {
			channel_configuration_list massless_configuration;
			channel_configuration_list massive_configuration;

			for (const sidis::partonic_channel::channel_type type : channel_types(pto)) {
				massless_configuration[type] = massless;
				massive_configuration[type] = massive;
			}

			return create_channels(pto, process, flavors, massless_configuration, massive_configuration);
		}

		/// @brief Builds the list of quark-to-quark scattering channels
		/// @param[in] process The scattering process, used to determine the sign of the gauge boson
		/// @param[in] parton_config Parton configuration, which determines the allowed flavors and the masses
		/// @param[in] massless A map of channel configurations for each SIDIS channel type, used in massless cases (~ no charm fragmentation).
		/// If `std::nullopt` is given as the value, that channel is not included.
		/// @param[in] massive A map of channel configurations for each SIDIS channel type, used in massive cases (~ massive charm fragmentation).
		/// If `std::nullopt` is given as the value, that channel is not included.
		static std::vector<sidis::partonic_channel> create_quark_to_quark_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &parton_config,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			// Output vector
			std::vector<sidis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const sidis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			// Allowed quark flavors, based on the sign of the gauge boson
			const std::vector<parton> &incoming_quarks = positive_W ? parton_config.allowed_lower_quarks : parton_config.allowed_upper_quarks;
			const std::vector<parton> &outgoing_quarks = positive_W ? parton_config.allowed_upper_quarks : parton_config.allowed_lower_quarks;
			const std::vector<parton> all_quark_flavors = collections::join(parton_config.allowed_lower_quarks, parton_config.allowed_upper_quarks);
		
			// Allowed antiquark flavors, based on the sign of the gauge boson
			const std::vector<parton> &incoming_antiquarks = positive_W ? parton_config.allowed_upper_antiquarks : parton_config.allowed_lower_antiquarks;
			const std::vector<parton> &outgoing_antiquarks = positive_W ? parton_config.allowed_lower_antiquarks : parton_config.allowed_upper_antiquarks;
			const std::vector<parton> all_antiquark_flavors = collections::join(parton_config.allowed_lower_antiquarks, parton_config.allowed_upper_antiquarks);

			// Charm and anticharm partons
			const parton c { .flavor = flavor::c, .mass = parton_config.mass(flavor::c) };
			const parton cbar { .flavor = flavor::cbar, .mass = parton_config.mass(flavor::cbar) };

			// Helper function for building quark-to-quark channels when
			// the final state contains exactly one quark (and potentially gluons)
			const auto quark_to_quark_channel = [&](
				const sidis::partonic_channel::channel_type type, 
				const parton &initial, 
				const parton &fragmenting, 
				const channel_configuration &config
			) {
				return sidis::partonic_channel {
					pto,
					type,
					initial,
					fragmenting,
					{},
					config.scheme,
					flavor::abs(fragmenting.flavor) == flavor::c ? std::vector{ fragmenting } : std::vector{ c, cbar }
				};
			};

			// Helper function for building quark-to-quark channels when
			// the final state contains exactly three quarks
			const auto quark_to_3quark_channel = [&](
				const sidis::partonic_channel::channel_type type, 
				const parton &initial, 
				const parton &fragmenting, 
				const parton &nonfragmenting1, 
				const parton &nonfragmenting2, 
				const channel_configuration &config
			) {
				return sidis::partonic_channel {
					pto,
					type,
					initial, fragmenting,
					{ nonfragmenting1, nonfragmenting2 },
					config.scheme,
					flavor::abs(fragmenting.flavor) == flavor::c ? std::vector{ fragmenting } : std::vector{ c, cbar }
				};
			};

			// Quark-to-quark channels
			for (const auto &[initial, fragmenting] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				if (const auto config = configuration(beta_alpha, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_quark_channel(beta_alpha, initial, fragmenting, *config));
				}
			}
			// Quark-to-quark flavor-conserving channels
			for (const auto &[initial, parton1, parton2] : std::views::cartesian_product(all_quark_flavors, incoming_quarks, outgoing_quarks)) {
				const parton fragmenting = initial;
				if (const auto config = configuration(qq_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qq_Fcon, initial, fragmenting, parton2, parton::conjugate_parton(parton1), *config));
				}
			}
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const parton fragmenting = initial;
				if (const auto config = configuration(alpha_alpha_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(alpha_alpha_Fcon, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Quark-to-antiquark channels
			for (const auto &[initial, fragmenting_conjugate] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(betabar_alpha, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(betabar_alpha, initial, fragmenting, fragmenting_conjugate, fragmenting_conjugate, *config));
				}
			}
			// Quark-to-antiquark flavor-conserving channels
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const parton fragmenting = parton::conjugate_parton(initial);
				if (const auto config = configuration(alphabar_alpha_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(alphabar_alpha_Fcon, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Quark-to-quark channels (arbitrary fragmenting quark)
			for (const auto &[initial, other, fragmenting] : std::views::cartesian_product(incoming_quarks, outgoing_quarks, all_quark_flavors)) {
				if (const auto config = configuration(qprime_alpha, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qprime_alpha, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Quark-to-quark channels (arbitrary initial quark)
			for (const auto &[initial, other_conjugate, fragmenting] : std::views::cartesian_product(all_quark_flavors, incoming_quarks, outgoing_quarks)) {
				const parton other = parton::conjugate_parton(other_conjugate);
				if (const auto config = configuration(betaprime_q, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(betaprime_q, initial, fragmenting, other, initial, *config));
				}
			}
			// Quark-to-antiquark channels (arbitrary fragmenting quark)
			for (const auto &[initial, other, fragmenting_conjugate] : std::views::cartesian_product(incoming_quarks, outgoing_quarks, all_quark_flavors)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(qbarprime_alpha, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qbarprime_alpha, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Quark-to-antiquark channels (arbitrary initial quark)
			for (const auto &[initial, other, fragmenting_conjugate] : std::views::cartesian_product(all_quark_flavors, outgoing_quarks, incoming_quarks)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(alphabarprime_q, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(alphabarprime_q, initial, fragmenting, other, initial, *config));
				}
			}

			// Antiquark-to-antiquark channels
			for (const auto &[initial, fragmenting] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				if (const auto config = configuration(alphabar_betabar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_quark_channel(alphabar_betabar, initial, fragmenting, *config));
				}
			}
			// Antiquark-to-antiquark flavor-conserving channels
			for (const auto &[initial, parton1, parton2] : std::views::cartesian_product(all_antiquark_flavors, incoming_antiquarks, outgoing_antiquarks)) {
				const parton fragmenting = initial;
				if (const auto config = configuration(qbar_qbar_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qbar_qbar_Fcon, initial, fragmenting, parton2, parton::conjugate_parton(parton1), *config));
				}
			}
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				const parton fragmenting = initial;
				if (const auto config = configuration(betabar_betabar_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(betabar_betabar_Fcon, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Antiquark-to-quark channels
			for (const auto &[initial, fragmenting_conjugate] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(alpha_betabar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(alpha_betabar, initial, fragmenting, fragmenting_conjugate, fragmenting_conjugate, *config));
				}
			}
			// Antiquark-to-quark flavor-conserving channels
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				const parton fragmenting = parton::conjugate_parton(initial);
				if (const auto config = configuration(beta_betabar_Fcon, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(beta_betabar_Fcon, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Antiquark-to-antiquark channels (arbitrary fragmenting antiquark)
			for (const auto &[initial, other, fragmenting] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks, all_antiquark_flavors)) {
				if (const auto config = configuration(qbarprime_betabar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qbarprime_betabar, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Antiquark-to-antiquark channels (arbitrary initial antiquark)
			for (const auto &[initial, other_conjugate, fragmenting] : std::views::cartesian_product(all_antiquark_flavors, incoming_antiquarks, outgoing_antiquarks)) {
				const parton other = parton::conjugate_parton(other_conjugate);
				if (const auto config = configuration(alphabarprime_qbar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(alphabarprime_qbar, initial, fragmenting, other, initial, *config));
				}
			}
			// Antiquark-to-quark channels (arbitrary fragmenting quark)
			for (const auto &[initial, other, fragmenting_conjugate] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks, all_antiquark_flavors)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(qprime_betabar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(qprime_betabar, initial, fragmenting, other, parton::conjugate_parton(fragmenting), *config));
				}
			}
			// Antiquark-to-quark channels (arbitrary initial antiquark)
			for (const auto &[initial, other, fragmenting_conjugate] : std::views::cartesian_product(all_antiquark_flavors, outgoing_antiquarks, incoming_antiquarks)) {
				const parton fragmenting = parton::conjugate_parton(fragmenting_conjugate);
				if (const auto config = configuration(betaprime_qbar, fragmenting.mass > 0.0)) {
					channels.push_back(quark_to_3quark_channel(betaprime_qbar, initial, fragmenting, other, initial, *config));
				}
			}

			return channels;
		}

		static std::vector<sidis::partonic_channel> create_quark_to_gluon_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			std::vector<sidis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const sidis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			const std::vector<parton> &incoming_quarks = positive_W ? flavors.allowed_lower_quarks : flavors.allowed_upper_quarks;
			const std::vector<parton> &outgoing_quarks = positive_W ? flavors.allowed_upper_quarks : flavors.allowed_lower_quarks;

			const std::vector<parton> &incoming_antiquarks = positive_W ? flavors.allowed_upper_antiquarks : flavors.allowed_lower_antiquarks;
			const std::vector<parton> &outgoing_antiquarks = positive_W ? flavors.allowed_lower_antiquarks : flavors.allowed_upper_antiquarks;

			const parton c { .flavor = flavor::c, .mass = flavors.mass(flavor::c) };
			const parton cbar { .flavor = flavor::cbar, .mass = flavors.mass(flavor::cbar) };

			const auto quark_to_gluon_channel = [&](
				const sidis::partonic_channel::channel_type type, 
				const parton &initial, 
				const parton &nonfragmenting_quark, 
				const channel_configuration &config
			) {
				return sidis::partonic_channel {
					pto,
					type,
					initial,
					parton { .flavor = flavor::g },
					{ nonfragmenting_quark },
					config.scheme,
					{ c, cbar }
				};
			};

			// Quark-to-gluon channels
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const bool massive_q = false;
				if (const auto config = configuration(g_alpha, massive_q)) {
					channels.push_back(quark_to_gluon_channel(g_alpha, initial, other, *config));
				}
			}
			// Antiquark-to-gluon channels
			for (const auto &[initial, other] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				const bool massive_q = false;
				if (const auto config = configuration(g_betabar, massive_q)) {
					channels.push_back(quark_to_gluon_channel(g_betabar, initial, other, *config));
				}
			}

			return channels;
		}
		static std::vector<sidis::partonic_channel> create_gluon_to_quark_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			std::vector<sidis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const sidis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			const std::vector<parton> &incoming_quarks = positive_W ? flavors.allowed_lower_quarks : flavors.allowed_upper_quarks;
			const std::vector<parton> &outgoing_quarks = positive_W ? flavors.allowed_upper_quarks : flavors.allowed_lower_quarks;

			const std::vector<parton> &incoming_antiquarks = positive_W ? flavors.allowed_upper_antiquarks : flavors.allowed_lower_antiquarks;
			const std::vector<parton> &outgoing_antiquarks = positive_W ? flavors.allowed_lower_antiquarks : flavors.allowed_upper_antiquarks;

			const parton c { .flavor = flavor::c, .mass = flavors.mass(flavor::c) };
			const parton cbar { .flavor = flavor::cbar, .mass = flavors.mass(flavor::cbar) };

			const auto gluon_to_quark_channel = [&](
				const sidis::partonic_channel::channel_type type, 
				const parton &nonfragmenting_quark, 
				const parton &fragmenting, 
				const channel_configuration &config
			) {
				return sidis::partonic_channel {
					pto,
					type,
					parton { .flavor = flavor::g },
					fragmenting,
					{ nonfragmenting_quark },
					config.scheme,
					flavor::abs(fragmenting.flavor) == flavor::c ? std::vector{ fragmenting } : std::vector{ c, cbar }
				};
			};

			// Gluon-to-quark channels
			for (const auto &[other_conjugate, fragmenting] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const parton other = parton::conjugate_parton(other_conjugate);
				if (const auto config = configuration(beta_g, fragmenting.mass > 0.0 || other.mass > 0.0)) {
					channels.push_back(gluon_to_quark_channel(beta_g, other, fragmenting, *config));
				}
			}
			// Gluon-to-antiquark channels
			for (const auto &[other_conjugate, fragmenting] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				const parton other = parton::conjugate_parton(other_conjugate);
				if (const auto config = configuration(alphabar_g, fragmenting.mass > 0.0 || other.mass > 0.0)) {
					channels.push_back(gluon_to_quark_channel(alphabar_g, other, fragmenting, *config));
				}
			}

			return channels;
		}
		static std::vector<sidis::partonic_channel> create_gluon_to_gluon_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			std::vector<sidis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const sidis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			const std::vector<parton> &incoming_quarks = positive_W ? flavors.allowed_lower_quarks : flavors.allowed_upper_quarks;
			const std::vector<parton> &outgoing_quarks = positive_W ? flavors.allowed_upper_quarks : flavors.allowed_lower_quarks;

			const parton c { .flavor = flavor::c, .mass = flavors.mass(flavor::c) };
			const parton cbar { .flavor = flavor::cbar, .mass = flavors.mass(flavor::cbar) };

			const auto gluon_to_gluon_channel = [&](const parton &nonfragmenting_quark1, const parton &nonfragmenting_quark2, const channel_configuration &config) {
				return sidis::partonic_channel {
					pto,
					gg,
					parton { .flavor = flavor::g },
					parton { .flavor = flavor::g },
					{ nonfragmenting_quark1, nonfragmenting_quark2 },
					config.scheme,
					std::vector{ c, cbar }
				};
			};

			// Gluon-to-gluon channels
			for (const auto &[parton1_conjugate, parton2] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				const parton parton1 = parton::conjugate_parton(parton1_conjugate);
				if (const auto config = configuration(gg, parton1.mass > 0.0 || parton2.mass > 0.0)) {
					channels.push_back(gluon_to_gluon_channel(parton1, parton2, *config));
				}
			}

			return channels;
		}
	};
}

}

#endif