#ifndef DM_SIDIS_COEFFICIENTS_HW_HELPER_IMPL_H
#define DM_SIDIS_COEFFICIENTS_HW_HELPER_IMPL_H

#include <cmath>										// std::log

#include "dimuon/common/splitting.hpp"

#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"

namespace dm {

namespace sidis::coefficients::_impl::helper::qq_nlo {
	inline vectorization::distribution_product unsubtracted(const sidis::evaluation_parameters &p) noexcept {
		const double xi = p.xi;
		const double zeta = p.zeta;

		const double m = p.channel.fragmentation_mass;
		const double m2 = math::power<2>(m);

		const double lambda = p.Q2 / (p.Q2 + m2);
		const double K_A = (1.0 - lambda) * std::log(1.0 - lambda) / lambda;

		const double zeta_p = ((1.0 - lambda) * xi + (1.0 - xi) * zeta) / (1.0 - lambda * xi);

		const double theta = static_cast<double>(zeta_p > p.z);

		const double double_delta_term = 0.5 * ((1.0 + 3.0 * lambda) * K_A + 1.0) / lambda + 4.0 + math::power<2>(std::numbers::pi) / 3.0;
		const vectorization::distribution_product double_delta = -p.pdf * p.ff * double_delta_term / (1.0 - p.pdf_momentum_fraction);
	
		const vectorization::distribution_product delta_splitting_1 = ((1.0 + math::power<2>(xi)) / (1.0 - xi)) * (p.pdf_hat - p.pdf);
		const vectorization::distribution_product delta_splitting_2 = p.pdf * (2.0 * std::log(1.0 - p.pdf_momentum_fraction) + 0.5 * p.pdf_momentum_fraction * (2.0 + p.pdf_momentum_fraction)) / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product delta_splitting = p.ff * (delta_splitting_1 + delta_splitting_2) * std::log((p.Q2 + m2) / p.factorization_scale);

		const double delta_regular_term = (1.0 - xi) * (1.0 + std::log((1.0 - xi) / xi)) - 2.0 * xi * std::log(xi) / (1.0 - xi);
		const vectorization::distribution_product delta_regular = delta_regular_term * p.pdf_hat * p.ff;

		const vectorization::distribution_product delta_distribution_term_1 = 2.0 * (std::log(1.0 - xi) / (1.0 - xi)) * (xi * p.pdf_hat - p.pdf);
		const vectorization::distribution_product delta_distribution_term_2 = p.pdf * math::power<2>(std::log(1.0 - p.pdf_momentum_fraction)) / (1.0 - p.pdf_momentum_fraction);
		const vectorization::distribution_product delta_distribution = (delta_distribution_term_1 + delta_distribution_term_2) * p.ff;

		const vectorization::distribution_product delta = delta_splitting + delta_regular + delta_distribution;

		const vectorization::distribution_product zeta_distribution = ((1.0 - xi) / (1.0 - zeta)) * p.pdf_hat * (p.ff_hat_p * theta - p.ff);

		const vectorization::distribution_product xi_distribution_1 = ((1.0 - xi) / math::power<2>(1.0 - lambda * xi)) * (p.pdf_hat * p.ff_hat_p * theta - p.pdf * p.ff);
		const vectorization::distribution_product xi_distribution_2 = p.pdf * p.ff * (((1.0 - lambda) * lambda * p.pdf_momentum_fraction) / (1.0 - lambda * p.pdf_momentum_fraction) + std::log(1.0 - lambda * p.pdf_momentum_fraction)) / (math::power<2>(lambda) * (1.0 - p.pdf_momentum_fraction));
		const vectorization::distribution_product xi_distribution = (1.0 - zeta) * (xi_distribution_1 + xi_distribution_2);

		const vectorization::distribution_product double_distribution = 2.0 * ((xi * p.pdf_hat * (zeta * p.ff_hat_p * theta - p.ff) / (1.0 - zeta) + p.pdf * p.ff) / (1.0 - xi) - p.pdf * p.ff * std::log(1.0 - p.pdf_momentum_fraction) / (1.0 - p.pdf_momentum_fraction));

		const vectorization::distribution_product regular_term = (1.0 - xi) * 2.0 * xi * zeta * p.pdf_hat * p.ff_hat_p * theta / (1.0 - lambda * xi);

		return double_delta + delta + zeta_distribution + xi_distribution + double_distribution + regular_term;;
	}
	
	inline vectorization::distribution_product subtraction_term(const sidis::evaluation_parameters &p) noexcept {
		if (p.zeta < p.z) { return 0.0; }

		const vectorization::distribution_product distribution_term_1 = -((1.0 + math::power<2>(p.zeta)) * (1.0 + 2.0 * std::log(1.0 - p.zeta))) * (p.ff_hat - p.ff) / (1.0 - p.zeta);
		const vectorization::distribution_product distribution_term_2 = (2.0 * p.z - std::log(1.0 - p.z) * (-1.0 + p.z * (2.0 + p.z) + 2.0 * std::log(1.0 - p.z))) * p.ff / (1.0 - p.z);
		const vectorization::distribution_product distribution_term = (distribution_term_1 + distribution_term_2) * p.pdf;

		const vectorization::distribution_product splitting_1 = ((1.0 + math::power<2>(p.zeta)) / (1.0 - p.zeta)) * (p.ff_hat - p.ff);
		const vectorization::distribution_product splitting_2 = p.ff * (2.0 * std::log(1.0 - p.z) + 0.5 * p.z * (2.0 + p.z)) / (1.0 - p.z);
		const vectorization::distribution_product splitting = p.pdf * (splitting_1 + splitting_2) * (std::log(p.fragmentation_scale) - 2.0 * std::log(p.channel.fragmentation_mass));

		return (distribution_term + splitting) / (1.0 - p.pdf_momentum_fraction);;
	}
}

namespace sidis::coefficients::_impl::helper::qg_nlo {
	inline vectorization::distribution_product lq_hq_common_term(const sidis::evaluation_parameters &p) noexcept {
		const double m = p.channel.fragmentation_mass;
		const double m2 = math::power<2>(m);

		const double lambda = p.Q2 / (p.Q2 + m2);

		const double zeta_min = std::max(p.z, ((1.0 - lambda) * p.xi) / (1.0 - lambda * p.xi));

		const auto delta = [&](const double xi) {
			const double splitting = splitting::Pqg(xi);
			const double log1 = std::log((p.Q2 + math::power<2>(p.channel.fragmentation_mass)) / p.factorization_scale);
			const double log2 = std::log((1.0 - xi) / xi);
			const double finite = xi * (1.0 - xi);

			return splitting * (log1 + log2) + finite;
		};
		const vectorization::distribution_product delta_term = p.pdf_hat * p.ff * delta(p.xi) / (1.0 - zeta_min);
		
		const vectorization::distribution_product distribution_term = p.pdf_hat * (((p.ff_hat - p.ff) / (1.0 - p.zeta) + p.ff_hat / p.zeta) + p.ff * std::log(1.0 - zeta_min) / (1.0 - zeta_min)) * splitting::Pqg(p.xi);

		return delta_term + distribution_term;
	}
}

}

#endif