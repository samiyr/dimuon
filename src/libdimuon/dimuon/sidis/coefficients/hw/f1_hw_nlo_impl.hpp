#ifndef DM_SIDIS_COEFFICIENTS_HW_F1_HW_NLO_IMPL_H
#define DM_SIDIS_COEFFICIENTS_HW_F1_HW_NLO_IMPL_H

#include "dimuon/common/constants.hpp"							// constants::T_R, constants::C_F
#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"
#include "dimuon/sidis/coefficients/hw/f2_hw_nlo_impl.hpp"
#include "dimuon/sidis/coefficients/hw/helper_impl.hpp"

namespace dm {

namespace sidis::coefficients::_impl::f1_hw::qq_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::qq_nlo::lq_lq(p);
		const vectorization::distribution_product F1_term = constants::C_F * p.pdf_hat * p.ff_hat * (4.0 * p.xi * p.zeta);
		return (F2_value / 2.0 - F1_term / p.z) / p.pdf_momentum_fraction;
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product unsubtracted_term = sidis::coefficients::_impl::helper::qq_nlo::unsubtracted(p);
		const vectorization::distribution_product subtraction_term = p.channel.include_subtraction_term(p.fragmentation_scale) ? sidis::coefficients::_impl::helper::qq_nlo::subtraction_term(p) : 0.0;

		return constants::C_F * (unsubtracted_term - subtraction_term) / (p.pdf_momentum_fraction * p.z);
	}
}

namespace sidis::coefficients::_impl::f1_hw::qg_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::qg_nlo::lq_lq(p);
		const vectorization::distribution_product F1_term = constants::T_R * p.pdf_hat * p.ff_hat * (8.0 * p.xi * (1.0 - p.xi));
		return (F2_value / 2.0 - F1_term / p.z) / p.pdf_momentum_fraction;
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const double m = p.channel.fragmentation_mass;
		const double m2 = math::power<2>(m);

		const double lambda = p.Q2 / (p.Q2 + m2);

		const double zeta_min = std::max(p.z, ((1.0 - lambda) * p.xi) / (1.0 - lambda * p.xi));

		if (p.zeta < zeta_min) { return 0.0; }

		const auto hg = [lambda](const double xi, const double zeta) {
			const double term1 = -math::power<2>(xi) * (1.0 - lambda) * (1.0 - 2.0 * lambda) / math::power<2>(zeta);
			const double term2 = 2.0 * xi * (1.0 - lambda) * (1.0 - 2.0 * lambda * xi) / zeta;
			const double term3 = 2.0 * xi * lambda * (1.0 - lambda * xi) - 1.0;

			return term1 + term2 + term3;
		};
		
		const vectorization::distribution_product common_term = sidis::coefficients::_impl::helper::qg_nlo::lq_hq_common_term(p);
		const vectorization::distribution_product hg_term = p.pdf_hat * p.ff_hat * hg(p.xi, p.zeta);

		return (common_term + hg_term) / (p.pdf_momentum_fraction * p.z);
	}
}

namespace sidis::coefficients::_impl::f1_hw::gq_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::gq_nlo::lq_lq(p);
		const vectorization::distribution_product F1_term = constants::C_F * p.pdf_hat * p.ff_hat * (4.0 * p.xi * (1.0 - p.zeta));
		return (F2_value / 2.0 - F1_term / p.z) / p.pdf_momentum_fraction;
	}	
}

}

#endif