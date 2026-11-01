#ifndef DM_DIS_PARTON_VALUES_H
#define DM_DIS_PARTON_VALUES_H

#include <array>
#include <cstddef>
#include <ranges>
#include <vector>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/flavor.hpp"

#include "dimuon/coupling/is_coupling.hpp"

#include "dimuon/pdf/is_pdf_interface.hpp"

#include "dimuon/dis/coefficient_channel.hpp"
#include "dimuon/dis/partonic_channel.hpp"
#include "dimuon/dis/parameters.hpp"
#include "dimuon/dis/vectorization.hpp"

namespace dm {

namespace dis {
	/// @brief PDF values at one integration point, shared by all partonic channels and structure functions
	/// @details Multiple channels need the same values, so each distinct value is computed once per integration point.
	/// The values are computed for all flavors at once.
	/// @tparam pdf_interface_t Type of PDF object. Must satisfy the `is_scalar_pdf_interface` concept.
	/// @tparam coupling_t Type of coupling object (e.g. `ckm_matrix`). Must satisfy the `is_coupling` concept.
	template <is_scalar_pdf_interface pdf_interface_t, is_coupling coupling_t>
	class parton_values {
	public:
		/// @brief Computes the PDF values at the given integration point
		/// @param[in] xi The convolution variable
		/// @param[in] params_ Parameters used for the evaluation of the DIS structure functions, which must outlive the values
		explicit parton_values(const double xi, const dis::parameters<pdf_interface_t, coupling_t> &params_) :
			params{params_},
			pdf_hat_values{
				params_.momentum_fractions
				| std::views::transform([&](const double momentum_fraction) {
					return compute_pdf_values(params_.pdf, momentum_fraction / xi, params_.factorization_scale);
				})
				| std::ranges::to<std::vector>()
			} {}

		/// @brief Returns the PDF values of \p channel at `momentum_fraction` and `momentum_fraction / xi` (`pdf` and `pdf_hat`)
		/// @param[in] channel The partonic channel (must be one of the channels used for the parameters)
		/// @param[in] coefficient_channel The coefficient channel of \p channel, which gives the PDF momentum fraction
		std::array<double, dis::vectorization::pdf_value_count> pdf_values(const dis::partonic_channel &channel, const dis::coefficient_channel &coefficient_channel) const {
			const std::size_t index = params.momentum_fraction_index(coefficient_channel.pdf_momentum_fraction(params.kinematics.x, params.kinematics.Q2));
			const std::size_t flavor_index = flavor::to_index(channel.initial_line.flavor);
			return {params.pdf_values_at_momentum_fractions[index][flavor_index], pdf_hat_values[index][flavor_index]};
		}

	private:
		/// The parameters of the integral
		const dis::parameters<pdf_interface_t, coupling_t> &params;
		/// The PDF values of all flavors at `momentum_fraction / xi`, for each momentum fraction in the order of `params.momentum_fractions`
		/// (the values at the momentum fractions themselves don't depend on the integration point, and are in `params`)
		const std::vector<std::array<double, constants::total_flavors>> pdf_hat_values;
	};
}

}

#endif
