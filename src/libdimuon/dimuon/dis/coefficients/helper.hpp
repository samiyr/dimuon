#ifndef DM_DIS_COEFFICIENTS_HELPER_H
#define DM_DIS_COEFFICIENTS_HELPER_H

#include <cmath>
#include <numbers>

#include "dimuon/common/splitting.hpp"
#include "dimuon/dis/coefficients/evaluation_parameters.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

namespace dis::coefficients::_impl::helper::gluon {
	/// The part of the massive gluon coefficient functions that is common to F1, F2, and F3. The sign of the log of `L` is \p sign: +1 for F1 and F2, and -1 for F3.
	inline double lq_hq_common_term(const dis::evaluation_parameters &p, const double sign) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);
		const double L = std::log((1.0 - lambda * p.xi) / ((1.0 - lambda) * p.xi));

		const double Pqg = splitting::Pqg(p.xi);

		const double splitting = Pqg * (sign * L + std::log((p.Q2 + m2) / p.factorization_scale));
		const double C0 = Pqg * (2.0 * std::log(1.0 - p.xi) - std::log(1.0 - lambda * p.xi) - std::log(p.xi));

		return splitting + C0;
	}

	/// The subtraction term of the massive gluon coefficient functions, the collinear log of the massive quark
	inline vectorization::distribution_product subtraction_term(const dis::evaluation_parameters &p) noexcept {
		const double m = p.channel.final_state_mass;
		const double log = std::log(p.factorization_scale / math::power<2>(m));
		const double splitting = splitting::Pqg(p.xi);

		return log * splitting * p.pdf_hat;
	}
}

namespace dis::coefficients::_impl::helper::quark {
	/// The part of the massive quark coefficient functions that is common to F1, F2, and F3
	inline vectorization::distribution_product lq_hq_common_term(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);
		const double K_A = (1.0 - lambda) * std::log(1.0 - lambda) / lambda;

		const double log = std::log((p.Q2 + m2) / p.factorization_scale);

		const vectorization::distribution_product splitting_term = ((1.0 + math::power<2>(p.xi)) / (1.0 - p.xi)) * (p.pdf_hat - p.pdf)
			+ p.pdf * (0.5 * p.pdf_momentum_fraction * (2.0 + p.pdf_momentum_fraction) + 2.0 * std::log(1.0 - p.pdf_momentum_fraction)) / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product splitting = log * splitting_term;

		const vectorization::distribution_product hq_delta = -p.pdf * (4.0 + 1.0 / (2.0 * lambda) + math::power<2>(std::numbers::pi) / 3.0 + K_A * (1.0 + 3.0 * lambda) / (2.0 * lambda))
			/ (1.0 - p.pdf_momentum_fraction);

		const vectorization::distribution_product hq_finite = -p.pdf_hat * (1.0 + math::power<2>(p.xi)) * std::log(p.xi) / (1.0 - p.xi);

		const vectorization::distribution_product hq_distribution_1 = ((2.0 * std::log(1.0 - p.xi) - std::log(1.0 - lambda * p.xi)) / (1.0 - p.xi))
			* ((1.0 + math::power<2>(p.xi)) * p.pdf_hat - 2.0 * p.pdf);
		const vectorization::distribution_product hq_distribution_2 = 2.0 * p.pdf * (
			std::log(1.0 - p.pdf_momentum_fraction) * std::log((1.0 - p.pdf_momentum_fraction) / (1.0 - lambda))
			- math::polylogT<2>(-lambda / (1.0 - lambda))
			+ math::polylogT<2>(-lambda * (1.0 - p.pdf_momentum_fraction) / (1.0 - lambda))
		) / (1.0 - p.pdf_momentum_fraction);

		const vectorization::distribution_product hq_distribution = hq_distribution_1 + hq_distribution_2;

		return splitting + hq_delta + hq_finite + hq_distribution;
	}

	/// The term `B3` of the massive quark coefficient functions, which is the same in F1, F2, and F3
	inline vectorization::distribution_product b3_term(const dis::evaluation_parameters &p) noexcept {
		const double m2 = math::power<2>(p.channel.final_state_mass);
		const double lambda = p.Q2 / (p.Q2 + m2);

		return 0.5 * (
			((1.0 - p.xi) / math::power<2>(1.0 - lambda * p.xi)) * (p.pdf_hat - p.pdf)
			+ p.pdf * (((1.0 - lambda) * lambda * p.pdf_momentum_fraction) / (1.0 - lambda * p.pdf_momentum_fraction) + std::log(1.0 - lambda * p.pdf_momentum_fraction))
				/ (math::power<2>(lambda) * (1.0 - p.pdf_momentum_fraction))
		);
	}
}

}

#endif
