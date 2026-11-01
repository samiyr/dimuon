#ifndef DM_SIDIS_COEFFICIENTS_EVALUATION_PARAMETERS_H
#define DM_SIDIS_COEFFICIENTS_EVALUATION_PARAMETERS_H

#include "dimuon/sidis/coefficients/evaluation_cache_gc.hpp"
#include "dimuon/sidis/coefficient_channel.hpp"
#include "dimuon/sidis/vectorization.hpp"

namespace dm {

namespace sidis {
	/// Struct of values needed for the evaluation of SIDIS coefficient functions.
	/// All values are for a specific kinematical point and integration point.
	/// The PDF and FF values are vectors with one lane per product of a PDF value and an FF value (see `sidis::vectorization::distribution_product`).
	/// The evaluator sets them to unit values (see `sidis::vectorization::make_masking_vector`), so that the coefficient functions give the coefficient of each product.
	struct evaluation_parameters {
		/// The properties of the SIDIS channel that determine its coefficient functions (see `coefficient_channel`). The coefficient functions are evaluated once
		/// for all the channels with the same coefficient channel, so they get only the coefficient channel, and not the channel with its flavors.
		const coefficient_channel &channel;

		/// Values of the gc coefficient functions that do not depend on the channel, computed once per integration point
		const evaluation_cache_gc &c;

		/// The PDF-side convolution variable
		double xi;
		/// PDF momentum fraction, either Bjorken-x or the slow-rescaling variable chi
		double pdf_momentum_fraction;

		/// Value (x*f(x)) of the appropriate PDF (based on the current channel) evaluated at `pdf_momentum_fraction`
		vectorization::distribution_product pdf;
		/// Value (x*f(x)) of the appropriate PDF (based on the current channel) evaluated at `pdf_momentum_fraction / xi`
		vectorization::distribution_product pdf_hat;

		/// The FF-side convolution variable
		double zeta;
		/// The FF-side momentum fraction
		double z;

		/// Value (z*D(z)) of the appropriate FF (based on the current channel) evaluated at `z`
		vectorization::distribution_product ff;
		/// Value (z*D(z)) of the appropriate FF (based on the current channel) evaluated at `z / zeta`
		vectorization::distribution_product ff_hat;
		/// Value (z*D(z)) of the appropriate FF (based on the current channel) evaluated at `z / zeta_p`,
		/// where `zeta_p` is related to `zeta` (see `sidis::parton_values::zeta_p`)
		vectorization::distribution_product ff_hat_p;

		/// The factorization scale at the current kinematical point
		double factorization_scale;
		/// The fragmentation scale at the current kinematical point
		double fragmentation_scale;

		/// The virtuality of the current kinematical point
		double Q2;
	};
}

}

#endif