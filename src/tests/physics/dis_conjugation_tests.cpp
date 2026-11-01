#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iterator>
#include <ranges>
#include <string>
#include <utility>

#include <gtest/gtest.h>

#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/perturbative_quantity.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/common/trf_kinematics.hpp>
#include <dimuon/dis/physical_channel.hpp>
#include <dimuon/dis/result.hpp>

#include "support/assertions.hpp"
#include "support/dis_setup.hpp"
#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum scattering_process::process_type;
	constexpr std::array orders{lo, nlo};

	/// The physical channel with the conjugate flavor of \p channel
	dis::physical_channel conjugate(dis::physical_channel channel) {
		channel.initial_flavor = flavor::conjugate_flavor(channel.initial_flavor);
		return channel;
	}

	std::string describe(const dis::physical_channel &channel) {
		return std::format("channel {}", std::to_underlying(channel.initial_flavor));
	}

	std::pair<dis::result, dis::result> conjugate_cross_sections() {
		const double charm_mass = 1.3;
		const auto antineutrino_dis = make_dis(antineutrino, mass_scheme_t::slow_rescaling, mass_scheme_t::sacot_subtracted, charm_mass, toy_pdf_interface());
		const auto neutrino_dis = make_dis(neutrino, mass_scheme_t::slow_rescaling, mass_scheme_t::sacot_subtracted, charm_mass, toy_pdf_interface(true));
		const trf_kinematics kinematics = trf_kinematics::Q2_s(0.05, 10.0, 300.0, antineutrino_dis.process.target.mass, antineutrino_dis.process.projectile.mass);
		return {antineutrino_dis.differential_cross_section_xy(kinematics), neutrino_dis.differential_cross_section_xy(kinematics)};
	}
}

TEST(dis_symmetry, charge_conjugation) {
	// Antineutrino scattering is the charge conjugate of neutrino scattering off the conjugate target
	const auto [antineutrino, conjugate_neutrino] = conjugate_cross_sections();
	for (const perturbative_order pto : orders) {
		EXPECT_REL_NEAR(antineutrino.value.correction(pto), conjugate_neutrino.value.correction(pto), 1e-14) << "order " << std::to_underlying(pto);
	}
}

TEST(dis_symmetry, channel_conjugation) {
	// Each physical channel is equal to the conjugate channel of the conjugate process
	const auto [antineutrino, conjugate_neutrino] = conjugate_cross_sections();
	ASSERT_EQ(antineutrino.physical_channels.size(), conjugate_neutrino.physical_channels.size());
	
	for (const auto &[channel, value] : std::views::zip(antineutrino.physical_channels, antineutrino.channel_decomposition)) {
		const auto found = conjugate_neutrino.physical_channels.find(conjugate(channel));
		ASSERT_NE(found, conjugate_neutrino.physical_channels.end()) << "no conjugate of the " << describe(channel);
		const perturbative_quantity &conjugate_value = conjugate_neutrino.channel_decomposition[static_cast<std::size_t>(std::distance(conjugate_neutrino.physical_channels.begin(), found))];
		for (const perturbative_order pto : orders) {
			EXPECT_NEAR(value.correction(pto), conjugate_value.correction(pto), 1e-10 * std::abs(antineutrino.value.correction(pto)))
				<< describe(channel) << ", order " << std::to_underlying(pto);
		}
	}
}
