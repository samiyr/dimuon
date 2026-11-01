#ifndef DM_SIDIS_COEFFICIENTS_HW_F2_HW_NLO_IMPL_H
#define DM_SIDIS_COEFFICIENTS_HW_F2_HW_NLO_IMPL_H

#include "dimuon/common/constants.hpp"							// constants::T_R, constants::C_F
#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"
#include "dimuon/sidis/coefficients/hw/helper_impl.hpp"
#include "dimuon/sidis/coefficients/hw/scale_impl.hpp"

namespace dm {

namespace sidis::coefficients::_impl::f2_hw::qq_nlo {
	inline vectorization::distribution_product unsubtracted(const sidis::evaluation_parameters &p) noexcept {
		const double xi = p.xi;
		const double zeta = p.zeta;

		const vectorization::distribution_product F1_value = sidis::coefficients::_impl::helper::qq_nlo::unsubtracted(p);

		const double lambda = p.Q2 / (p.Q2 + math::power<2>(p.channel.fragmentation_mass));
		const double K_A = (1.0 - lambda) * std::log(1.0 - lambda) / lambda;

		const double zeta_p = ((1.0 - lambda) * xi + (1.0 - xi) * zeta) / (1.0 - lambda * xi);

		const double theta = static_cast<double>(zeta_p > p.z);

		const vectorization::distribution_product double_delta = K_A * p.pdf * p.ff / (1.0 - p.pdf_momentum_fraction);

		const double regular_term_1 = -2.0 * zeta * xi * (1.0 - 3.0 * lambda);
		const double regular_term_2 = -2.0 * (1.0 - lambda);
		const vectorization::distribution_product regular_term = p.pdf_hat * p.ff_hat_p * theta * (regular_term_1 + regular_term_2) * (1.0 - xi) / (1.0 - lambda * xi);

		const vectorization::distribution_product F2_term = double_delta + regular_term;

		return F1_value + F2_term;
	}

	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const vectorization::distribution_product term1_delta = p.pdf * p.ff * (math::power<(2)>((std::log(1.0 - p.pdf_momentum_fraction) + std::log(1.0 - p.z))) - 8.0);

		const vectorization::distribution_product term2_xi = p.pdf_hat * p.ff * ((1.0 - p.xi) * (1.0 + std::log((1.0 - p.xi) / p.xi) + std::log(1.0 - p.z)) - (2.0 * p.xi * std::log(p.xi)) / (1.0 - p.xi));
		const vectorization::distribution_product term3_xi = (2.0 * (std::log(1.0 - p.xi) + std::log(1.0 - p.z))) * ((p.xi * p.pdf_hat * p.ff - p.pdf * p.ff) / (1.0 - p.xi));

		const vectorization::distribution_product term4_zeta = p.pdf * p.ff_hat * ((1.0 - p.zeta) * (1.0 + std::log(p.zeta * (1.0 - p.zeta)) + std::log(1.0 - p.pdf_momentum_fraction)) + (2.0 * p.zeta * std::log(p.zeta)) / (1.0 - p.zeta));
		const vectorization::distribution_product term5_zeta = (2.0 * (std::log(1.0 - p.zeta) + std::log(1.0 - p.pdf_momentum_fraction))) * ((p.zeta * p.pdf * p.ff_hat - p.pdf * p.ff) / (1.0 - p.zeta));

		const vectorization::distribution_product term6_xi_zeta = (1.0 - p.xi) / (1.0 - p.zeta) * (p.pdf_hat * p.ff_hat - p.pdf_hat * p.ff);
		const vectorization::distribution_product term7_xi_zeta = (1.0 - p.zeta) / (1.0 - p.xi) * (p.pdf_hat * p.ff_hat - p.pdf * p.ff_hat);
		const vectorization::distribution_product term8_xi_zeta = 2.0 * (p.xi * p.zeta * p.pdf_hat * p.ff_hat - p.xi * p.pdf_hat * p.ff - p.zeta * p.pdf * p.ff_hat + p.pdf * p.ff) / ((1.0 - p.xi) * (1.0 - p.zeta));
		const vectorization::distribution_product term9_xi_zeta = 6.0 * p.xi * p.zeta * p.pdf_hat * p.ff_hat;

		const vectorization::distribution_product result = term1_delta / ((1.0 - p.pdf_momentum_fraction) * (1.0 - p.z)) + (term2_xi + term3_xi) / (1.0 - p.z) + (term4_zeta + term5_zeta) / (1.0 - p.pdf_momentum_fraction) + (term6_xi_zeta + term7_xi_zeta + term8_xi_zeta + term9_xi_zeta);
		const vectorization::distribution_product coefficient_contribution = 2.0 * constants::C_F * result / p.z;

		const vectorization::distribution_product scale_contribution = sidis::coefficients::_impl::scale::quark_to_quark(p);
		
		return coefficient_contribution + scale_contribution;
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		const vectorization::distribution_product unsubtracted_term = sidis::coefficients::_impl::f2_hw::qq_nlo::unsubtracted(p);
		vectorization::distribution_product subtraction_term = p.channel.include_subtraction_term(p.fragmentation_scale) ? sidis::coefficients::_impl::helper::qq_nlo::subtraction_term(p) : 0.0;

		return 2.0 * constants::C_F * (unsubtracted_term - subtraction_term) / p.z;
	}
}

namespace sidis::coefficients::_impl::f2_hw::qg_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const vectorization::distribution_product term1_xi = p.pdf_hat * p.ff * (1.0 - (math::power<2>(p.xi) + math::power<2>(1.0 - p.xi)) * (1.0 - std::log((1.0 - p.xi) / p.xi)));
		const vectorization::distribution_product term2_xi = p.pdf_hat * p.ff * (std::log(1.0 - p.z) * (1.0 - 2.0 * p.xi * (1.0 - p.xi)));

		const vectorization::distribution_product term3_xi_zeta = p.pdf_hat * p.ff_hat * (12.0 * p.xi * (1.0 - p.xi) + (1.0 - p.zeta) / p.zeta);
		const vectorization::distribution_product term4_xi_zeta = (p.pdf_hat * p.ff_hat * (p.zeta - 2.0 * p.xi * (1.0 - p.xi) / p.zeta) - p.pdf_hat * p.ff * (1.0 - 2.0 * p.xi * (1.0 - p.xi))) / (1.0 - p.zeta);

		const vectorization::distribution_product result = (term1_xi + term2_xi) / (1.0 - p.z) + (term3_xi_zeta + term4_xi_zeta);
		const vectorization::distribution_product coefficient_contribution = 2.0 * constants::T_R * result / p.z;

		const vectorization::distribution_product scale_contribution = sidis::coefficients::_impl::scale::gluon_to_quark(p);
					
		return coefficient_contribution + scale_contribution;
	}

	inline vectorization::distribution_product lq_hq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const double m = p.channel.fragmentation_mass;
		const double m2 = math::power<2>(m);

		const double lambda = p.Q2 / (p.Q2 + m2);

		const double zeta_min = std::max(p.z, ((1.0 - lambda) * p.xi) / (1.0 - lambda * p.xi));

		if (p.zeta < zeta_min) { return 0.0; }

		const auto hg = [lambda](const double xi, const double zeta) {
			const double term1 = math::power<2>(xi) * (1.0 - lambda) * (1.0 - 6.0 * lambda + 6.0 * math::power<2>(lambda)) / math::power<2>(zeta);
			const double term2 = 6.0 * lambda * xi * (1.0 - lambda) * (1.0 - 2.0 * lambda * xi) / zeta;
			const double term3 = lambda * (6.0 * lambda * xi * (1.0 - lambda * xi) - 1.0);

			return term1 + term2 + term3;
		};

		const vectorization::distribution_product common_term = sidis::coefficients::_impl::helper::qg_nlo::lq_hq_common_term(p);		
		const vectorization::distribution_product hg_term = p.pdf_hat * p.ff_hat * hg(p.xi, p.zeta);

		return 2.0 * (common_term + hg_term) / p.z;
	}
}

namespace sidis::coefficients::_impl::f2_hw::gq_nlo {
	inline vectorization::distribution_product lq_lq(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const vectorization::distribution_product term1_zeta = p.pdf * p.ff_hat * (p.zeta + std::log(p.zeta * (1.0 - p.zeta)) * (1 + math::power<2>(1.0 - p.zeta)) / p.zeta);
		const vectorization::distribution_product term2_zeta = p.pdf * p.ff_hat * (std::log(1.0 - p.pdf_momentum_fraction) * (p.zeta + 2.0 * (1.0 - p.zeta) / p.zeta));

		const vectorization::distribution_product term3_xi_zeta = p.pdf_hat * p.ff_hat * (6.0 * p.xi * (1.0 - p.zeta) + (1.0 - p.xi) / p.zeta);
		const vectorization::distribution_product term4_xi_zeta = (p.pdf_hat * p.ff_hat * (p.zeta + 2.0 * p.xi * (1.0 - p.zeta) / p.zeta) - p.pdf * p.ff_hat * (p.zeta + 2.0 * (1.0 - p.zeta) / p.zeta)) / (1.0 - p.xi);

		const vectorization::distribution_product result = (term1_zeta + term2_zeta) / (1.0 - p.pdf_momentum_fraction) + (term3_xi_zeta + term4_xi_zeta);
		const vectorization::distribution_product coefficient_contribution = 2.0 * constants::C_F * result / p.z;

		const vectorization::distribution_product scale_contribution = sidis::coefficients::_impl::scale::quark_to_gluon(p);
		
		return coefficient_contribution + scale_contribution;
	}	
}

}

#endif