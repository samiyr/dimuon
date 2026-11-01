#include <algorithm>
#include <array>
#include <ranges>
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/dis/channel_factory.hpp>
#include <dimuon/dis/partonic_channel.hpp>

using namespace dm;

namespace {
	using enum perturbative_order;
	using enum scattering_process::process_type;

	const std::vector<flavor::type> flavors{flavor::u, flavor::d, flavor::s, flavor::c, flavor::b};

	scattering_process make_process(const scattering_process::process_type type) {
		return scattering_process{type, constants::particles::proton, constants::particles::neutrino};
	}

	/// All channels of the order \p pto in \p process, with different mass schemes for the massless and the massive (charm and bottom) channels
	std::vector<dis::partonic_channel> make_channels(const perturbative_order pto, const scattering_process::process_type process) {
		std::array<double, constants::total_flavors> masses{};
		for (const flavor::type f : {flavor::c, flavor::cbar}) { masses[flavor::to_index(f)] = 1.3; }
		for (const flavor::type f : {flavor::b, flavor::bbar}) { masses[flavor::to_index(f)] = 4.75; }
		return dis::channel_factory::create_all_channels(
			pto, make_process(process), parton_configuration{flavors, masses},
			channel_configuration{mass_scheme_t::slow_rescaling}, channel_configuration{mass_scheme_t::sacot_subtracted}
		);
	}
}

TEST(dis_channels, charge_conservation) {
	// Check that every channel conserves charge
	for (const auto process : {neutrino, antineutrino}) {
		for (const perturbative_order pto : {lo, nlo}) {
			for (const dis::partonic_channel &channel : make_channels(pto, process)) {
				EXPECT_TRUE(channel.check_charge_conservation(make_process(process))) << channel;
			}
		}
	}
}

TEST(dis_channels, coupling) {
	// Check that the coupling is always non-zero
	const ckm_matrix V = ckm_matrix::global_fit();
	for (const auto process : {neutrino, antineutrino}) {
		for (const perturbative_order pto : {lo, nlo}) {
			for (const dis::partonic_channel &channel : make_channels(pto, process)) {
				const auto [a, b] = channel.get_coupling_quarks();
				EXPECT_GT(V(a.flavor, b.flavor), 0.0) << channel;
			}
		}
	}
}

TEST(dis_channels, channel_conjugation) {
	// Check that the antineutrino channels are the charge conjugates of the neutrino channels, up to the order of the final-state quarks
	for (const perturbative_order pto : {lo, nlo}) {
		std::vector<dis::partonic_channel> conjugated_neutrino = make_channels(pto, neutrino);
		for (dis::partonic_channel &channel : conjugated_neutrino) {
			channel.type = dis::partonic_channel::conjugate_type(channel.type);
			channel.initial_line = parton::conjugate_parton(channel.initial_line);
			for (parton &p : channel.final_lines) { p = parton::conjugate_parton(p); }
			std::ranges::sort(channel.final_lines, {}, &parton::flavor);
		}
		std::vector<dis::partonic_channel> antineutrino_channels = make_channels(pto, antineutrino);
		for (dis::partonic_channel &channel : antineutrino_channels) { std::ranges::sort(channel.final_lines, {}, &parton::flavor); }
		EXPECT_TRUE(std::ranges::is_permutation(conjugated_neutrino, antineutrino_channels)) << "order " << std::to_underlying(pto);
	}
}