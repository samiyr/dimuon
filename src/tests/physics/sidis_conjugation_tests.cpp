#include <cmath>
#include <format>
#include <ranges>
#include <utility>

#include <gtest/gtest.h>

#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/perturbative_quantity.hpp>
#include <dimuon/sidis/physical_channel.hpp>
#include <dimuon/sidis/result.hpp>

#include "support/assertions.hpp"
#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum scattering_process::process_type;
	constexpr std::array orders{lo, nlo, nnlo};

	/// The physical channel with the conjugate flavors of \p channel
	sidis::physical_channel conjugate(sidis::physical_channel channel) {
		channel.initial_flavor = flavor::conjugate_flavor(channel.initial_flavor);
		channel.fragmentation_flavor = flavor::conjugate_flavor(channel.fragmentation_flavor);
		return channel;
	}

	std::string describe(const sidis::physical_channel &channel) {
		return std::format("channel {} --> {}", std::to_underlying(channel.initial_flavor), std::to_underlying(channel.fragmentation_flavor));
	}

	std::pair<sidis::result, sidis::result> conjugate_cross_sections() {
		return {
			compute_sidis(sidis_run{.process = antineutrino}, toy_pdf_interface(), toy_ff_interface()),
			compute_sidis(sidis_run{.process = neutrino}, toy_pdf_interface(true), toy_ff_interface(true))
		};
	}
}

TEST(sidis_symmetry, charge_conjugation) {
	// Antineutrino scattering is the charge conjugate of neutrino scattering off the conjugate target into the conjugate hadrons
	const auto [antineutrino, conjugate_neutrino] = conjugate_cross_sections();
	for (const perturbative_order pto : orders) {
		EXPECT_REL_NEAR(antineutrino.value.correction(pto), conjugate_neutrino.value.correction(pto), 1e-12) << "order " << std::to_underlying(pto);
	}
}

TEST(sidis_symmetry, channel_conjugation) {
	// Each physical channel is equal to the conjugate channel of the conjugate process
	const auto [antineutrino, conjugate_neutrino] = conjugate_cross_sections();
	ASSERT_EQ(antineutrino.physical_channels.size(), conjugate_neutrino.physical_channels.size());
	for (const auto &[index, channel] : std::views::enumerate(antineutrino.physical_channels)) {
		const auto found = conjugate_neutrino.physical_channels.find(conjugate(channel));
		ASSERT_NE(found, conjugate_neutrino.physical_channels.end()) << "no conjugate of the " << describe(channel);
		const auto conjugate_index = std::distance(conjugate_neutrino.physical_channels.begin(), found);
		for (const perturbative_order pto : orders) {
			EXPECT_NEAR(
				(antineutrino.channel_decomposition[0, 0, index].correction(pto)), (conjugate_neutrino.channel_decomposition[0, 0, conjugate_index].correction(pto)),
				1e-12 * std::abs(antineutrino.value.correction(pto))
			) << describe(channel) << ", order " << std::to_underlying(pto);
		}
	}
}
