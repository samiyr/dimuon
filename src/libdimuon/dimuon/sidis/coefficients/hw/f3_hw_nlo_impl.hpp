#ifndef DM_SIDIS_COEFFICIENTS_HW_F3_HW_NLO_IMPL_H
#define DM_SIDIS_COEFFICIENTS_HW_F3_HW_NLO_IMPL_H

#include "dimuon/common/constants.hpp"							// constants::T_R, constants::C_F
#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"
#include "dimuon/sidis/coefficients/hw/f2_hw_nlo_impl.hpp"
#include "dimuon/sidis/coefficients/hw/helper_impl.hpp"

namespace dm {

namespace sidis::coefficients::_impl::f3_hw::qq_nlo {
	inline vectorization::distribution_product unsubtracted(const sidis::evaluation_parameters &p) noexcept {
		const double xi = p.xi;
		const double zeta = p.zeta;

		const vectorization::distribution_product F1_value = sidis::coefficients::_impl::helper::qq_nlo::unsubtracted(p);

		const double lambda = p.Q2 / (p.Q2 + math::power<2>(p.channel.fragmentation_mass));

		const double zeta_p = ((1.0 - lambda) * xi + (1.0 - xi) * zeta) / (1.0 - lambda * xi);

		const double theta = static_cast<double>(zeta_p > p.z);

		const double regular_term_1 = 2.0 * zeta * (1.0 - xi);
		const double regular_term_2 = -2.0 * (1.0 - lambda * xi);
		const vectorization::distribution_product regular_term = p.pdf_hat * p.ff_hat_p * theta * (regular_term_1 + regular_term_2) * (1.0 - xi) / (1.0 - lambda * xi);

		const vectorization::distribution_product F3_term = regular_term;

		return F1_value + F3_term;
	}

	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::qq_nlo::lq_lq(p);
		const vectorization::distribution_product F3_term = constants::C_F * p.pdf_hat * p.ff_hat * (6.0 * p.xi * p.zeta + 2.0 * (1.0 - p.xi - p.zeta));
		return F2_value / p.pdf_momentum_fraction - 2.0 * F3_term / (p.pdf_momentum_fraction * p.z);
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product unsubtracted_term = sidis::coefficients::_impl::f3_hw::qq_nlo::unsubtracted(p);
		const vectorization::distribution_product subtraction_term = p.channel.include_subtraction_term(p.fragmentation_scale) ? sidis::coefficients::_impl::helper::qq_nlo::subtraction_term(p) : 0.0;

		return 2.0 * constants::C_F * (unsubtracted_term - subtraction_term) / (p.pdf_momentum_fraction * p.z);
	}
}

namespace sidis::coefficients::_impl::f3_hw::qg_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::qg_nlo::lq_lq(p);
		const vectorization::distribution_product F3_term = constants::T_R * p.pdf_hat * p.ff_hat * (12.0 * p.xi * (1.0 - p.xi) + 2.0 * (1.0 - 2.0 * p.xi * (1.0 - p.xi)) / p.zeta - 2.0);
		return F2_value / p.pdf_momentum_fraction - 2.0 * F3_term / (p.pdf_momentum_fraction * p.z);
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const double m = p.channel.fragmentation_mass;
		const double m2 = math::power<2>(m);

		const double lambda = p.Q2 / (p.Q2 + m2);

		const double zeta_min = std::max(p.z, ((1.0 - lambda) * p.xi) / (1.0 - lambda * p.xi));

		if (p.zeta < zeta_min) { return 0.0; }

		const auto hg = [lambda](const double xi, const double zeta) {
			const double splitting = splitting::Pqg(xi);

			const double term1 = math::power<2>(xi) * (1.0 - lambda) * (1.0 - 2.0 * lambda) / math::power<2>(zeta);
			const double term2 = -2.0 * (splitting + xi * (1.0 - xi - lambda * xi) * (1.0 - lambda)) / zeta;

			return term1 + term2;
		};

		const vectorization::distribution_product common_term = sidis::coefficients::_impl::helper::qg_nlo::lq_hq_common_term(p);		
		const vectorization::distribution_product hg_term = p.pdf_hat * p.ff_hat * hg(p.xi, p.zeta);

		return 2.0 * (common_term + hg_term) / (p.pdf_momentum_fraction * p.z);
	}
}

namespace sidis::coefficients::_impl::f3_hw::gq_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		// The F2 coefficient function without the division of the modified schemes, which applies to F2 only
		const vectorization::distribution_product F2_value = sidis::coefficients::_impl::f2_hw::gq_nlo::lq_lq(p);
		const vectorization::distribution_product F3_term = constants::C_F * p.pdf_hat * p.ff_hat * (4.0 * p.xi * (1.0 - p.zeta) + 2.0 * (1.0 - p.xi) * p.zeta);
		return F2_value / p.pdf_momentum_fraction - 2.0 * F3_term / (p.pdf_momentum_fraction * p.z);
	}
}

}

#endif