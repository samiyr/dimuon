#ifndef DM_DIS_COEFFICIENTS_F1_H
#define DM_DIS_COEFFICIENTS_F1_H

#include <cmath>
#include <utility>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/dis/coefficients/evaluation_parameters.hpp"
#include "dimuon/dis/coefficients/f2.hpp"
#include "dimuon/dis/coefficients/helper.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

namespace dis::coefficients::_impl::f1::gluon_nlo {
	inline vectorization::distribution_product lq_lq(const dis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product F2_value = dis::coefficients::_impl::f2::gluon_nlo::lq_lq(p);
		const vectorization::distribution_product F1_term = p.pdf_hat * 4.0 * p.xi * (1.0 - p.xi);
		return (F2_value - 2.0 * F1_term) / (2.0 * p.pdf_momentum_fraction);
	}

	inline vectorization::distribution_product unsubtracted(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);
		const double L = std::log((1.0 - lambda * p.xi) / ((1.0 - lambda) * p.xi));

		const double C1 = 4.0 - 4.0 * (1.0 - lambda);
		const double C2 = ((1.0 - lambda) * p.xi) / (1.0 - lambda * p.xi) - 1.0;
		const double C3 = 2.0;
		const double C4 = -4.0;

		const double common = dis::coefficients::_impl::helper::gluon::lq_hq_common_term(p, 1.0);
		const double term = C1 * p.xi * (1.0 - p.xi) + C2 + (1.0 - lambda) * p.xi * L * (C3 + lambda * p.xi * C4);

		return p.pdf_hat * (common + term);
	}

	inline vectorization::distribution_product lq_hq(const dis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product unsubtracted_term = dis::coefficients::_impl::f1::gluon_nlo::unsubtracted(p);
		const vectorization::distribution_product subtraction_term = p.channel.include_subtraction_term(p.factorization_scale)
			? dis::coefficients::_impl::helper::gluon::subtraction_term(p)
			: 0.0;

		return (unsubtracted_term - subtraction_term) / p.pdf_momentum_fraction;
	}
}

namespace dis::coefficients::_impl::f1::quark_nlo {
	inline vectorization::distribution_product lq_lq(const dis::evaluation_parameters &p) noexcept {
		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = dis::coefficients::_impl::f2::quark_nlo::lq_lq(p);
		const vectorization::distribution_product F1_term = constants::C_F * p.pdf_hat * (2.0 * p.xi);
		return (F2_value - 2.0 * F1_term) / (2.0 * p.pdf_momentum_fraction);
	}

	inline vectorization::distribution_product lq_hq(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);

		const vectorization::distribution_product common = dis::coefficients::_impl::helper::quark::lq_hq_common_term(p);

		const vectorization::distribution_product B1 = ((1.0 - 4.0 * p.xi + math::power<2>(p.xi)) * p.pdf_hat + 2.0 * p.pdf) / (1.0 - p.xi)
			- 2.0 * p.pdf * std::log(1.0 - p.pdf_momentum_fraction) / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product B2 = p.pdf_hat * p.xi * (1.0 - p.xi) / (1.0 - lambda * p.xi);
		const vectorization::distribution_product B3 = dis::coefficients::_impl::helper::quark::b3_term(p);

		const vectorization::distribution_product term = B1 + B2 + B3;

		return constants::C_F * (common + term) / p.pdf_momentum_fraction;
	}
}

/// Namespace containing functions for the evaluation of the F1 structure function
namespace dis::coefficients::f1 {
	namespace lo {
		/// @brief Evaluates the LO quark F1 structure function
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark(const dis::evaluation_parameters &p) noexcept {
			// p.pdf gives x*f(x), so divide x out. Also, divide out the (at-LO useless) integration over xi.
			// No mass scheme choice at LO, since the slow-rescaling prescription is exact at LO
			return p.pdf / p.pdf_momentum_fraction / (1.0 - p.pdf_momentum_fraction);
		}

		/// @brief Evaluates the trivial LO gluon F1 structure function
		/// @param[in] p Evaluation parameter is ignored
		/// @return `0.0`
		/// @note This function exists only to keep the function structure as exists at higher orders
		inline vectorization::distribution_product gluon(const dis::evaluation_parameters &) noexcept { return 0.0; }
	}

	namespace nlo {
		/// @brief Evaluates the NLO quark F1 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark(const dis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return dis::coefficients::_impl::f1::quark_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return dis::coefficients::_impl::f1::quark_nlo::lq_hq(p);
			}
			std::unreachable();
		}

		/// @brief Evaluates the NLO gluon F1 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product gluon(const dis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return dis::coefficients::_impl::f1::gluon_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return dis::coefficients::_impl::f1::gluon_nlo::lq_hq(p);
			}
			std::unreachable();
		}
	}
}

}

#endif
