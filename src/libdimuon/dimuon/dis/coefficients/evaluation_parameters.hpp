#ifndef DM_DIS_COEFFICIENTS_EVALUATION_PARAMETERS_H
#define DM_DIS_COEFFICIENTS_EVALUATION_PARAMETERS_H

#include "dimuon/dis/coefficient_channel.hpp"
#include "dimuon/dis/vectorization.hpp"

namespace dm {

namespace dis {
	/// Struct of values needed for the evaluation of DIS coefficient functions.
	/// All values are for a specific kinematical point and integration point.
	/// The PDF values are vectors with one lane per PDF value (see `dis::vectorization::distribution_product`).
	/// The evaluator sets them to unit values (see `dis::vectorization::make_masking_vector`), so that the coefficient functions give the coefficient of each PDF value.
	struct evaluation_parameters {
		/// The properties of the DIS channel that determine its coefficient functions (see `coefficient_channel`). The coefficient functions are evaluated once
		/// for all the channels with the same coefficient channel, so they get only the coefficient channel, and not the channel with its flavors.
		const coefficient_channel &channel;

		/// The convolution variable
		double xi;
		/// PDF momentum fraction, either Bjorken-x or the slow-rescaling variable chi
		double pdf_momentum_fraction;

		/// Value (x*f(x)) of the appropriate PDF (based on the current channel) evaluated at `pdf_momentum_fraction`
		vectorization::distribution_product pdf;
		/// Value (x*f(x)) of the appropriate PDF (based on the current channel) evaluated at `pdf_momentum_fraction / xi`
		vectorization::distribution_product pdf_hat;

		/// The factorization scale (GeV^2) at the current kinematical point
		double factorization_scale;

		/// The virtuality (GeV^2) of the current kinematical point
		double Q2;
	};
}

}

#endif
