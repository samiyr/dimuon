#ifndef DM_DIS_COEFFICIENTS_F2_H
#define DM_DIS_COEFFICIENTS_F2_H

#include <cmath>
#include <numbers>
#include <utility>

#include "dimuon/common/constants.hpp"
#include "dimuon/common/mass_scheme.hpp"
#include "dimuon/dis/coefficients/evaluation_parameters.hpp"
#include "dimuon/dis/coefficients/helper.hpp"
#include "dimuon/dis/coefficients/scale.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

namespace dis::coefficients::_impl::f2::gluon_nlo {
	inline vectorization::distribution_product lq_lq(const dis::evaluation_parameters &p) noexcept {
		const double term1_xi = (math::power<2>(p.xi) + math::power<2>(1.0 - p.xi)) * std::log((1.0 - p.xi) / p.xi);
		const double term2_xi = -1.0 + 8.0 * p.xi * (1.0 - p.xi);

		const vectorization::distribution_product coefficient_contribution = 2.0 * p.pdf_hat * (term1_xi + term2_xi);

		const vectorization::distribution_product scale_contribution = dis::coefficients::_impl::scale::gluon(p);

		return coefficient_contribution + scale_contribution;
	}

	inline vectorization::distribution_product unsubtracted(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);
		const double L = std::log((1.0 - lambda * p.xi) / ((1.0 - lambda) * p.xi));

		const double C1 = 8.0 - 18.0 * (1.0 - lambda) + 12.0 * math::power<2>(1.0 - lambda);
		const double C2 = (1.0 - lambda) / (1.0 - lambda * p.xi) - 1.0;
		const double C3 = 6.0 * lambda;
		const double C4 = -12.0 * lambda;

		const double common = dis::coefficients::_impl::helper::gluon::lq_hq_common_term(p, 1.0);
		const double term = C1 * p.xi * (1.0 - p.xi) + C2 + (1.0 - lambda) * p.xi * L * (C3 + lambda * p.xi * C4);

		return p.pdf_hat * (common + term);
	}

	inline vectorization::distribution_product lq_hq(const dis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product unsubtracted_term = dis::coefficients::_impl::f2::gluon_nlo::unsubtracted(p);
		vectorization::distribution_product subtraction_term = p.channel.include_subtraction_term(p.factorization_scale)
			? dis::coefficients::_impl::helper::gluon::subtraction_term(p)
			: 0.0;

		if (p.channel.is_modified_mass_scheme()) {
			subtraction_term /= p.channel.chi_factor(p.Q2);
		}

		return 2.0 * (unsubtracted_term - subtraction_term);
	}
}

namespace dis::coefficients::_impl::f2::quark_nlo {
	inline vectorization::distribution_product lq_lq(const dis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product term1_delta = -p.pdf * (
			9.0 / 2.0 + math::power<2>(std::numbers::pi) / 3.0 - math::power<2>(std::log(1.0 - p.pdf_momentum_fraction)) + 1.5 * std::log(1.0 - p.pdf_momentum_fraction)
		);

		const vectorization::distribution_product term2_xi = std::log(1.0 - p.xi) * ((1.0 + math::power<2>(p.xi)) * p.pdf_hat - 2.0 * p.pdf) / (1.0 - p.xi);
		const vectorization::distribution_product term3_xi = (p.pdf_hat - p.pdf) / (1.0 - p.xi);
		const vectorization::distribution_product term4_xi = p.pdf_hat * (-(1.0 + math::power<2>(p.xi)) * std::log(p.xi) / (1.0 - p.xi) + 3.0 + 2.0 * p.xi);

		const vectorization::distribution_product result = term1_delta / (1.0 - p.pdf_momentum_fraction) + (term2_xi - 1.5 * term3_xi + term4_xi);
		const vectorization::distribution_product coefficient_contribution = 2.0 * constants::C_F * result;

		const vectorization::distribution_product scale_contribution = dis::coefficients::_impl::scale::quark(p);

		return coefficient_contribution + scale_contribution;
	}

	inline vectorization::distribution_product lq_hq(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);
		const double K_A = (1.0 - lambda) * std::log(1.0 - lambda) / lambda;

		const vectorization::distribution_product common = dis::coefficients::_impl::helper::quark::lq_hq_common_term(p);

		const vectorization::distribution_product A = K_A * p.pdf / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product B1 = ((2.0 - 2.0 * math::power<2>(p.xi) - 2.0 / p.xi) * p.pdf_hat + 2.0 * p.pdf) / (1.0 - p.xi)
			- 2.0 * p.pdf * std::log(1.0 - p.pdf_momentum_fraction) / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product B2 = p.pdf_hat * (2.0 / p.xi - 1.0 - p.xi) / (1.0 - lambda * p.xi);
		const vectorization::distribution_product B3 = dis::coefficients::_impl::helper::quark::b3_term(p);

		const vectorization::distribution_product term = A + B1 + B2 + B3;

		return 2.0 * constants::C_F * (common + term);
	}
}

/// Namespace containing functions for the evaluation of the F2 structure function
namespace dis::coefficients::f2 {
	namespace lo {
		/// @brief Evaluates the LO quark F2 structure function
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark(const dis::evaluation_parameters &p) noexcept {
			// p.pdf gives x*f(x); the (at-LO useless) integration over xi is divided out
			vectorization::distribution_product value = 2.0 * p.pdf / (1.0 - p.pdf_momentum_fraction);
			if (p.channel.is_modified_mass_scheme()) {
				value /= p.channel.chi_factor(p.Q2);
			}
			return value;
		}

		/// @brief Evaluates the trivial LO gluon F2 structure function
		/// @param[in] p Evaluation parameter is ignored
		/// @return `0.0`
		/// @note This function exists only to keep the function structure as exists at higher orders
		inline vectorization::distribution_product gluon(const dis::evaluation_parameters &) noexcept { return 0.0; }
	}

	namespace nlo {
		/// @brief Evaluates the NLO quark F2 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product quark(const dis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
				return dis::coefficients::_impl::f2::quark_nlo::lq_lq(p);
			case slow_rescaling_modified:
				return dis::coefficients::_impl::f2::quark_nlo::lq_lq(p) / p.channel.chi_factor(p.Q2);
			case sacot_unsubtracted:
			case sacot_subtracted:
				return dis::coefficients::_impl::f2::quark_nlo::lq_hq(p);
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return dis::coefficients::_impl::f2::quark_nlo::lq_hq(p) / p.channel.chi_factor(p.Q2);
			}
			std::unreachable();
		}

		/// @brief Evaluates the NLO gluon F2 structure function in the requested mass scheme
		/// @param[in] p Evaluation parameters
		inline vectorization::distribution_product gluon(const dis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return dis::coefficients::_impl::f2::gluon_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return dis::coefficients::_impl::f2::gluon_nlo::lq_hq(p);
			}
			std::unreachable();
		}
	}
}

}

#endif
