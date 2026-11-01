#ifndef DM_EXE_CHANNELS_H
#define DM_EXE_CHANNELS_H

#include <vector>
#include <format>
#include <stdexcept>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/flavor.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/dis/partonic_channel.hpp>
#include <dimuon/sidis/partonic_channel.hpp>
#include <dimuon/dis/channel_factory.hpp>
#include <dimuon/sidis/channel_factory.hpp>

#include "maps.hpp"
#include <argparser/arg_parser.hpp>

/// Reads the channels from the arguments.
/// A channel type takes the mass scheme of its own argument, e.g. --nlo-ba-m for the massless alpha --> beta SIDIS channels at NLO
/// (and --nlo-ba-M for the massive ones), or --nlo-q-m for the massless quark DIS channels, or else that of the argument for all channel types
/// of the order, e.g. --nlo-all-m.
namespace channels {
	namespace {
		std::optional<dm::channel_configuration> get_mass_scheme_arg(const arg_parser &args, const std::string arg) {
			const std::optional<std::string> arg_opt = args.get_argument(arg);
			if (!arg_opt) return std::nullopt;
			
			const std::string &arg_value = *arg_opt;
			if (!maps::mass_schemes.contains(arg_value)) {
				throw std::runtime_error(std::format("Invalid mass scheme '{}' of argument --{}", arg_value, arg));
			}

			return std::make_optional<dm::channel_configuration>(dm::channel_configuration{ .scheme = maps::mass_schemes.at(arg_value) });
		}

		/// The configuration of the channel type \p type at the order \p order, see the description of the namespace
		std::optional<dm::channel_configuration> get_channel_configuration(const arg_parser &args, const std::string order, const std::string type, const std::string mass_type) {
			return get_mass_scheme_arg(args, std::format("{}-{}-{}", order, type, mass_type))
				.or_else([&] { return get_mass_scheme_arg(args, std::format("{}-all-{}", order, mass_type)); });
		}
	}

	inline std::vector<dm::dis::partonic_channel> read_dis_channels(
		const std::string order, 
		const dm::parton_configuration &flavors, 
		const dm::scattering_process &process, 
		const arg_parser &args
	) {
		using enum dm::dis::partonic_channel::channel_type;

		const dm::perturbative_order pto = maps::perturbative_orders.at(order);

		// The names of the channel types in the arguments
		const std::unordered_map<dm::dis::partonic_channel::channel_type, std::string> type_names = {
			{beta_alpha,		"q"},
			{alphabar_betabar,	"qbar"},
			{beta_g,			"g"},
			{alphabar_g,		"gbar"}
		};

		// The configurations of the channel types of the order (the gluon channels start at NLO)
		const auto create_channel_map = [&](const std::string &mass_type) {
			dm::dis::channel_factory::channel_configuration_list channel_map;
			for (const dm::dis::partonic_channel::channel_type type : dm::dis::channel_factory::channel_types(pto)) {
				channel_map[type] = get_channel_configuration(args, order, type_names.at(type), mass_type);
			}
			return channel_map;
		};

		return dm::dis::channel_factory::create_channels(
			pto, process, flavors,
			create_channel_map("m"),
			create_channel_map("M")
		);
	}

	inline std::vector<dm::sidis::partonic_channel> read_sidis_channels(
		const std::string order, 
		const dm::parton_configuration &flavors, 
		const dm::scattering_process &process, 
		const arg_parser &args
	) {
		using enum dm::sidis::partonic_channel::channel_type;

		const dm::perturbative_order pto = maps::perturbative_orders.at(order);

		// The names of the channel types in the arguments
		const std::unordered_map<dm::sidis::partonic_channel::channel_type, std::string> type_names = {
			{beta_alpha,			"ba"},
			{qq_Fcon,				"qqFcon"},
			{alpha_alpha_Fcon,		"aaFcon"},
			{betabar_alpha,			"bba"},
			{alphabar_alpha_Fcon,	"aba"},
			{qprime_alpha,			"qpa"},
			{betaprime_q,			"bpq"},
			{qbarprime_alpha,		"qbpa"},
			{alphabarprime_q,		"abpq"},
				
			{alphabar_betabar,		"abbb"},
			{qbar_qbar_Fcon,		"qbqbFcon"},
			{betabar_betabar_Fcon,	"bbbbFcon"},
			{alpha_betabar,			"abb"},
			{beta_betabar_Fcon,		"bbbFcon"},
			{qprime_betabar,		"qpbb"},
			{betaprime_qbar,		"bpqb"},
			{qbarprime_betabar,		"qbpbb"},
			{alphabarprime_qbar,	"abpqb"},
				
			{g_alpha,				"ga"},
			{g_betabar,				"gbb"},
				
			{beta_g,				"bg"},
			{alphabar_g,			"abg"},

			{gg,					"gg"}
		};

		// The configurations of the channel types of the order
		const auto create_channel_map = [&](const std::string &mass_type) {
			dm::sidis::channel_factory::channel_configuration_list channel_map;
			for (const dm::sidis::partonic_channel::channel_type type : dm::sidis::channel_factory::channel_types(pto)) {
				channel_map[type] = get_channel_configuration(args, order, type_names.at(type), mass_type);
			}
			return channel_map;
		};

		return dm::sidis::channel_factory::create_channels(
			pto, process, flavors,
			create_channel_map("m"),
			create_channel_map("M")
		);
	}

	constexpr inline std::vector<dm::sidis::partonic_channel> filter_channels(
		const std::vector<dm::sidis::partonic_channel> &channels, 
		const dm::flavor::flavor_vector &initial_lines, 
		const dm::flavor::flavor_vector &fragmentation_lines
	) {
		return channels | std::views::filter([&](const dm::sidis::partonic_channel &channel) {
			return 	std::ranges::contains(initial_lines, channel.initial_line.flavor) 
				&&  std::ranges::contains(fragmentation_lines, channel.fragmentation_line.flavor);
		}) | std::ranges::to<std::vector>();
	}
}

#endif