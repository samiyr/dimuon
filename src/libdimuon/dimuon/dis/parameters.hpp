#ifndef DM_DIS_PARAMETERS_H
#define DM_DIS_PARAMETERS_H

#include <array>
#include <cstddef>
#include <flat_set>
#include <ranges>
#include <utility>
#include <vector>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/coupling/is_coupling.hpp"
#include "dimuon/pdf/is_pdf_interface.hpp"
#include "dimuon/dis/coefficient_channel.hpp"
#include "dimuon/dis/partonic_channel.hpp"
#include "dimuon/utility/collections.hpp"

namespace dm {

namespace dis {
	/// @brief Returns the values of the PDF at \p momentum_fraction for all flavors, or zero if \p momentum_fraction is above one
	/// (i.e. where PDFs vanish). This gives the lower bound of the convolution integral, `momentum_fraction / xi <= 1` for `xi >= momentum_fraction`.
	/// @param[in] pdf PDF
	/// @param[in] momentum_fraction Momentum fraction
	/// @param[in] scale Momentum scale (GeV^2)
	std::array<double, constants::total_flavors> compute_pdf_values(const is_scalar_pdf_interface auto &pdf, const double momentum_fraction, const double scale) {
		return momentum_fraction <= 1.0 ? pdf.xfxQ2(momentum_fraction, scale) : std::array<double, constants::total_flavors>{};
	}

	/// Container for holding parameters related to the entire DIS structure function calculation
	template <is_scalar_pdf_interface pdf_interface_t, is_coupling coupling_t>
	struct parameters {
		/// PDF
		pdf_interface_t pdf;
		/// Coupling
		coupling_t coupling;
		/// Scattering process
		scattering_process process;
		/// Kinematical point
		trf_kinematics kinematics;
		/// Factorization scale (GeV^2) at the kinematical point `kinematics`
		double factorization_scale;
		/// Set of unique PDF momentum fractions (i.e. x and all possible values of the slow-rescaling variable across all channels)
		std::flat_set<double> momentum_fractions;
		/// The PDF values of all flavors at each of `momentum_fractions`, in the same order (see `compute_pdf_values`)
		std::vector<std::array<double, constants::total_flavors>> pdf_values_at_momentum_fractions;

		explicit parameters(
			pdf_interface_t pdf_,
			coupling_t coupling_,
			scattering_process process_,
			trf_kinematics kinematics_,
			const std::vector<dis::partonic_channel> &channels,
			double factorization_scale_
		) :
		pdf{std::move(pdf_)},
		coupling{std::move(coupling_)},
		process{std::move(process_)},
		kinematics{std::move(kinematics_)},
		factorization_scale{factorization_scale_},
		momentum_fractions{
			channels
			| std::views::transform([&](const dis::partonic_channel &channel) {
				return dis::coefficient_channel{channel}.pdf_momentum_fraction(kinematics.x, kinematics.Q2);
			})
			| std::ranges::to<std::flat_set<double>>()
		},
		pdf_values_at_momentum_fractions{
			momentum_fractions
			| std::views::transform([&](const double momentum_fraction) { return compute_pdf_values(pdf, momentum_fraction, factorization_scale); })
			| std::ranges::to<std::vector>()
		} {}

		/// @brief Returns the position of \p momentum_fraction in `momentum_fractions`
		/// @details It is assumed that \p momentum_fraction is contained in `momentum_fractions`
		std::size_t momentum_fraction_index(const double momentum_fraction) const noexcept {
			return collections::iterator_to_index(momentum_fractions, momentum_fractions.find(momentum_fraction));
		}
	};
}

}

#endif
