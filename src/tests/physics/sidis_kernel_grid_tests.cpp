#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/perturbative_quantity.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/interpolation/linear_interpolation_basis.hpp>
#include <dimuon/pdf/interface/functional_form_interface.hpp>
#include <dimuon/pdf/interface/vector_kernel_interface.hpp>
#include <dimuon/sidis/physical_channel.hpp>
#include <dimuon/sidis/result.hpp>

#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum scattering_process::process_type;
	constexpr std::array orders{lo, nlo, nnlo};

	std::string describe(const sidis::physical_channel &channel) {
		return std::format("channel {} --> {}", std::to_underlying(channel.initial_flavor), std::to_underlying(channel.fragmentation_flavor));
	}

	/// Piecewise-linear interpolation of \p value(k) on \p knots, zero outside them, written independently of `linear_interpolation_basis`
	double interpolate(const std::vector<double> &knots, const auto &value, const double x) {
		if (x < knots.front() || x > knots.back()) { return 0.0; }
		const std::size_t right = std::min(static_cast<std::size_t>(std::ranges::upper_bound(knots, x) - knots.begin()), knots.size() - 1);
		const double t = (x - knots[right - 1]) / (knots[right] - knots[right - 1]);
		return (1.0 - t) * value(right - 1) + t * value(right);
	}
}

TEST(sidis_kernel_grid, convolution) {
	// A kernel grid holds the cross section of each pair of interpolation kernels. For a PDF and an FF that are piecewise linear on the knots,
	// the grid times their values at the knots is their cross section, in each physical channel. With plain QMC the integrands are equal at every point.
	const std::vector<double> pdf_knots{0.08, 0.12, 0.2, 0.35, 0.5, 0.7, 0.9, 1.0};
	const std::vector<double> ff_knots{0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0};
	const functional_form_interface piecewise_linear_pdf{[&](const flavor::type f, const double x, const double) {
		return interpolate(pdf_knots, [&](const std::size_t k) { return toy_pdf(f, pdf_knots[k]); }, x);
	}};
	const functional_form_interface piecewise_linear_ff{[&](const flavor::type f, const double z, const double) {
		return interpolate(ff_knots, [&](const std::size_t k) { return toy_ff(f, ff_knots[k]); }, z);
	}};

	for (const auto process : {neutrino, antineutrino}) {
		SCOPED_TRACE(std::format("process {}", std::to_underlying(process)));
		const sidis::result direct = compute_sidis(sidis_run{.process = process}, piecewise_linear_pdf, piecewise_linear_ff);
		const sidis::result grid = compute_sidis(
			sidis_run{.process = process}, vector_kernel_interface{linear_interpolation_basis{pdf_knots}}, vector_kernel_interface{linear_interpolation_basis{ff_knots}}
		);
		ASSERT_EQ(grid.channel_decomposition.extent(0), pdf_knots.size());
		ASSERT_EQ(grid.channel_decomposition.extent(1), ff_knots.size());
		ASSERT_TRUE(std::ranges::equal(grid.physical_channels, direct.physical_channels));

		for (const perturbative_order pto : orders) {
			for (const auto &[c, channel] : std::views::enumerate(grid.physical_channels)) {
				double convolution = 0.0;
				for (std::size_t k = 0; k < pdf_knots.size(); ++k) {
					for (std::size_t l = 0; l < ff_knots.size(); ++l) {
						convolution += grid.channel_decomposition[k, l, c].correction(pto) * toy_pdf(channel.initial_flavor, pdf_knots[k]) * toy_ff(channel.fragmentation_flavor, ff_knots[l]);
					}
				}
				EXPECT_NEAR(convolution, (direct.channel_decomposition[0, 0, c].correction(pto)), 1e-10 * std::abs(direct.value.correction(pto)))
					<< describe(channel) << ", order " << std::to_underlying(pto);
			}
		}
	}
}
