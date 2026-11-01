#ifndef DM_DIS_COEFFICIENTS_SCALE_H
#define DM_DIS_COEFFICIENTS_SCALE_H

#include <cmath>

#include "dimuon/common/constants.hpp"
#include "dimuon/dis/coefficients/evaluation_parameters.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

/// The factorization-scale dependence of the massless coefficient functions
namespace dis::coefficients::_impl::scale {
	inline vectorization::distribution_product quark(const dis::evaluation_parameters &p) noexcept {
		const double factorization_scale_log = std::log(p.Q2 / p.factorization_scale);

		const vectorization::distribution_product term1_delta = factorization_scale_log * p.pdf * (2.0 * std::log(1.0 - p.pdf_momentum_fraction) + 1.5);
		const vectorization::distribution_product term2_xi = factorization_scale_log * (p.pdf_hat * (1.0 + math::power<2>(p.xi)) - 2.0 * p.pdf) / (1.0 - p.xi);

		const vectorization::distribution_product result = term1_delta / (1.0 - p.pdf_momentum_fraction) + term2_xi;
		return 2.0 * constants::C_F * result;
	}

	inline vectorization::distribution_product gluon(const dis::evaluation_parameters &p) noexcept {
		const double factorization_scale_log = std::log(p.Q2 / p.factorization_scale);

		return 2.0 * factorization_scale_log * p.pdf_hat * (math::power<2>(1.0 - p.xi) + math::power<2>(p.xi));
	}
}

}

#endif
