#ifndef DM_SIDIS_COEFFICIENTS_HW_SCALE_IMPL_H
#define DM_SIDIS_COEFFICIENTS_HW_SCALE_IMPL_H

#include "dimuon/common/constants.hpp"

#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"

namespace dm {

namespace sidis::coefficients::_impl::scale {
	inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
		const double factorization_scale_log = std::log(p.Q2 / p.factorization_scale);
		const double fragmentation_scale_log = std::log(p.Q2 / p.fragmentation_scale);

		const vectorization::distribution_product term1_delta = p.pdf * p.ff * ((2.0 * std::log(1.0 - p.z) + 1.5) * fragmentation_scale_log + (2.0 * std::log(1.0 - p.pdf_momentum_fraction) + 1.5) * factorization_scale_log);
		const vectorization::distribution_product term2_xi = factorization_scale_log * (((1.0 + math::power<2>(p.xi)) * p.pdf_hat * p.ff - 2.0 * p.pdf * p.ff) / (1.0 - p.xi));
		const vectorization::distribution_product term3_zeta = fragmentation_scale_log * (((1.0 + math::power<2>(p.zeta)) * p.pdf * p.ff_hat - 2.0 * p.pdf * p.ff) / (1.0 - p.zeta));

		const vectorization::distribution_product result = term1_delta / ((1.0 - p.pdf_momentum_fraction) * (1.0 - p.z)) + term2_xi / (1.0 - p.z) + term3_zeta / (1.0 - p.pdf_momentum_fraction);
		return 2.0 * constants::C_F * result / p.z;
	}

	inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &p) noexcept {
		const double fragmentation_scale_log = std::log(p.Q2 / p.fragmentation_scale);

		const vectorization::distribution_product result = fragmentation_scale_log * p.pdf * p.ff_hat * ((1.0 + math::power<2>(1.0 - p.zeta)) / p.zeta) / (1.0 - p.pdf_momentum_fraction);
		return 2.0 * constants::C_F * result / p.z;
	}

	inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &p) noexcept {
		const double factorization_scale_log = std::log(p.Q2 / p.factorization_scale);

		const vectorization::distribution_product result = factorization_scale_log * p.pdf_hat * p.ff * ((math::power<2>(p.xi) + math::power<2>(1.0 - p.xi))) / (1.0 - p.z);
		return 2.0 * constants::T_R * result / p.z;
	}
}

}

#endif