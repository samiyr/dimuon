#ifndef DM_DIS_CHANNEL_FACTORY_H
#define DM_DIS_CHANNEL_FACTORY_H

#include <optional>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include "dimuon/common/channel_configuration.hpp"
#include "dimuon/common/flavor.hpp"
#include "dimuon/common/parton.hpp"
#include "dimuon/common/parton_configuration.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/dis/partonic_channel.hpp"
#include "dimuon/utility/collections.hpp"

namespace dm {

namespace dis {
	/// Factory class for constructing DIS channels, with separation for massless (~ no charm) and massive (~ charm) channels
	struct channel_factory {
		using enum dis::partonic_channel::channel_type;

		/// The configurations of the channel types: a channel type that is missing, or whose configuration is `std::nullopt`, is not included
		using channel_configuration_list = std::unordered_map<dis::partonic_channel::channel_type, std::optional<channel_configuration>>;

		/// Returns the types of the channels at the perturbative order \p pto. There are no NNLO DIS coefficient functions, and so no NNLO channels.
		static std::vector<dis::partonic_channel::channel_type> channel_types(const perturbative_order pto) {
			switch (pto) {
				case perturbative_order::lo:
					return { beta_alpha, alphabar_betabar };
				case perturbative_order::nlo:
					return { beta_alpha, alphabar_betabar, beta_g, alphabar_g };
				case perturbative_order::nnlo:
					return {};
			}
			std::unreachable();
		}

		/// @brief Builds the list of channels of all types: the quark and antiquark channels, followed by the gluon channels
		/// @param[in] pto Perturbative order
		/// @param[in] process The scattering process, used to determine the sign of the gauge boson
		/// @param[in] flavors Parton configuration, which determines the allowed flavors and the masses
		/// @param[in] massless A map of channel configurations for each DIS channel type, used in massless cases (~ no charm).
		/// If `std::nullopt` is given as the value, or the channel type is missing, that channel is not included.
		/// @param[in] massive A map of channel configurations for each DIS channel type, used in massive cases (~ massive charm).
		/// If `std::nullopt` is given as the value, or the channel type is missing, that channel is not included.
		static std::vector<dis::partonic_channel> create_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			return collections::join(
				create_quark_channels(pto, process, flavors, massless, massive),
				create_gluon_channels(pto, process, flavors, massless, massive)
			);
		}

		/// Builds the list of all channels of the perturbative order \p pto (see `channel_types`), with the configuration \p massless
		/// for the massless channels and \p massive for the massive channels. If a configuration is `std::nullopt`, those channels are not included.
		static std::vector<dis::partonic_channel> create_all_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const std::optional<channel_configuration> &massless,
			const std::optional<channel_configuration> &massive
		) {
			channel_configuration_list massless_configuration;
			channel_configuration_list massive_configuration;

			for (const dis::partonic_channel::channel_type type : channel_types(pto)) {
				massless_configuration[type] = massless;
				massive_configuration[type] = massive;
			}

			return create_channels(pto, process, flavors, massless_configuration, massive_configuration);
		}

		/// @brief Builds the list of quark and antiquark channels. A channel is massive if its final-state quark is massive.
		/// @param[in] pto Perturbative order
		/// @param[in] process The scattering process, used to determine the sign of the gauge boson
		/// @param[in] flavors Parton configuration, which determines the allowed flavors and the masses
		/// @param[in] massless A map of channel configurations for each DIS channel type, used in massless cases (see `create_channels`)
		/// @param[in] massive A map of channel configurations for each DIS channel type, used in massive cases (see `create_channels`)
		static std::vector<dis::partonic_channel> create_quark_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			std::vector<dis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const dis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			// Allowed quark and antiquark flavors, based on the sign of the gauge boson
			const std::vector<parton> &incoming_quarks = positive_W ? flavors.allowed_lower_quarks : flavors.allowed_upper_quarks;
			const std::vector<parton> &outgoing_quarks = positive_W ? flavors.allowed_upper_quarks : flavors.allowed_lower_quarks;
			const std::vector<parton> &incoming_antiquarks = positive_W ? flavors.allowed_upper_antiquarks : flavors.allowed_lower_antiquarks;
			const std::vector<parton> &outgoing_antiquarks = positive_W ? flavors.allowed_lower_antiquarks : flavors.allowed_upper_antiquarks;

			// Adds the channel of the given type, with the massive configuration if the final-state quark is massive, and the massless one otherwise
			const auto add_channel = [&](const dis::partonic_channel::channel_type type, const parton &initial, const parton &final) {
				if (const auto config = configuration(type, final.mass > 0.0)) {
					channels.push_back(dis::partonic_channel{pto, type, initial, {final}, config->scheme});
				}
			};

			for (const auto &[initial, final] : std::views::cartesian_product(incoming_quarks, outgoing_quarks)) {
				add_channel(beta_alpha, initial, final);
			}
			for (const auto &[initial, final] : std::views::cartesian_product(incoming_antiquarks, outgoing_antiquarks)) {
				add_channel(alphabar_betabar, initial, final);
			}

			return channels;
		}

		/// @brief Builds the list of gluon channels: the gauge boson turns an (anti)quark of a gluon splitting into another one, which gives one channel for each
		/// final-state pair. For a W+ these are the channels of type `beta_g`, and for a W- the charge-conjugate channels of type `alphabar_g`.
		/// A channel is massive if either of its final-state quarks is massive.
		/// @param[in] pto Perturbative order
		/// @param[in] process The scattering process, used to determine the sign of the gauge boson
		/// @param[in] flavors Parton configuration, which determines the allowed flavors and the masses
		/// @param[in] massless A map of channel configurations for each DIS channel type, used in massless cases (see `create_channels`)
		/// @param[in] massive A map of channel configurations for each DIS channel type, used in massive cases (see `create_channels`)
		static std::vector<dis::partonic_channel> create_gluon_channels(
			const perturbative_order pto,
			const scattering_process &process,
			const parton_configuration &flavors,
			const channel_configuration_list &massless,
			const channel_configuration_list &massive
		) {
			std::vector<dis::partonic_channel> channels;

			// Returns the configuration of a channel of type `type`, from `massive` if the channel is massive and from `massless` otherwise,
			// or `std::nullopt` if the channel is not included (the type is missing, or its configuration is `std::nullopt`)
			const auto configuration = [&](const dis::partonic_channel::channel_type type, const bool is_massive) -> std::optional<channel_configuration> {
				const auto &configurations = is_massive ? massive : massless;
				const auto entry = configurations.find(type);
				return entry != configurations.end() ? entry->second : std::nullopt;
			};

			const bool positive_W = process.positive_W();

			// The (anti)quarks that the gauge boson turns into each other: a gluon splits into the conjugate of the incoming one and the outgoing one
			const dis::partonic_channel::channel_type type = positive_W ? beta_g : alphabar_g;
			const std::vector<parton> &incoming = positive_W ? flavors.allowed_lower_quarks : flavors.allowed_lower_antiquarks;
			const std::vector<parton> &outgoing = positive_W ? flavors.allowed_upper_quarks : flavors.allowed_upper_antiquarks;

			for (const auto &[other_conjugate, final] : std::views::cartesian_product(incoming, outgoing)) {
				if (const auto config = configuration(type, other_conjugate.mass > 0.0 || final.mass > 0.0)) {
					channels.push_back(dis::partonic_channel{
						pto, type, parton{ .flavor = flavor::g }, {final, parton::conjugate_parton(other_conjugate)}, config->scheme
					});
				}
			}

			return channels;
		}
	};
}

}

#endif
