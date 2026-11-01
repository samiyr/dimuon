#ifndef DM_SIDIS_COEFFICIENTS_F1_HW_H
#define DM_SIDIS_COEFFICIENTS_F1_HW_H

#include <utility>

#include "dimuon/sidis/coefficients/hw/f1_hw_nlo_impl.hpp"

namespace dm {

/// Namespace containing functions for the evaluation of the "handwritten" F1 structure function
namespace sidis::coefficients::f1_hw {
	namespace lo {
		/// @brief Evaluates the LO quark-to-quark F1 structure function
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			// p.pdf gives x*f(x) and p.ff gives z*D(z), so divide those out
			// Also, divide out the (at-LO useless) integration over xi
			// No mass scheme choice at LO, since the slow-rescaling prescription is exact at LO
			return p.pdf * p.ff / (p.pdf_momentum_fraction * p.z) / (1.0 - p.pdf_momentum_fraction);
		}

		/// @brief Evaluates the trivial LO quark-to-gluon F1 structure function
		/// @param[in] p Evaluation parameter is ignored
		/// @return `0.0`
		/// @note This function exists only to keep the function structure as exists at higher orders
		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }

		/// @brief Evaluates the trivial LO gluon-to-quark F1 structure function
		/// @param[in] p Evaluation parameter is ignored
		/// @return `0.0`
		/// @note This function exists only to keep the function structure as exists at higher orders
		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &) noexcept { return 0.0; }

		/// @brief Evaluates the trivial LO gluon-to-gluon F1 structure function
		/// @param[in] p Evaluation parameter is ignored
		/// @return `0.0`
		/// @note This function exists only to keep the function structure as exists at higher orders
		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}

	namespace nlo {
		/// @brief Evaluates the NLO quark-to-quark F1 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			// Dispatch to the appropriate implementation based on the mass scheme
			
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return sidis::coefficients::_impl::f1_hw::qq_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return sidis::coefficients::_impl::f1_hw::qq_nlo::lq_hq(p);
			}
			std::unreachable();
		}

		/// @brief Evaluates the NLO quark-to-gluon F1 structure function
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &p) noexcept {
			// No dispatch, in all schemes the quark-to-gluon contribution is computed with
			// massless quarks (and the slow-rescaling prescription)
			return sidis::coefficients::_impl::f1_hw::gq_nlo::lq_lq(p);
		}

		/// @brief Evaluates the NLO gluon-to-quark F1 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &p) noexcept {
			// Dispatch to the appropriate implementation based on the mass scheme

			using enum mass_scheme_t;

			// TODO check schemes

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return sidis::coefficients::_impl::f1_hw::qg_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return sidis::coefficients::_impl::f1_hw::qg_nlo::lq_hq(p);
			}
			std::unreachable();
		}

		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}
}

}

#endif