#ifndef DM_SIDIS_COEFFICIENTS_F3_HW_H
#define DM_SIDIS_COEFFICIENTS_F3_HW_H

#include <utility>

#include "dimuon/sidis/coefficients/hw/f3_hw_nlo_impl.hpp"

namespace dm {

namespace sidis::coefficients::f3_hw {
	namespace lo {
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return 2.0 * p.pdf * p.ff / (p.pdf_momentum_fraction * p.z) / (1.0 - p.pdf_momentum_fraction);
		}

		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &) noexcept { return 0.0; }
		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}

	namespace nlo {
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return sidis::coefficients::_impl::f3_hw::qq_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return sidis::coefficients::_impl::f3_hw::qq_nlo::lq_hq(p);
			}
			std::unreachable();
		}

		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &p) noexcept {
			return sidis::coefficients::_impl::f3_hw::gq_nlo::lq_lq(p);
		}

		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &p) noexcept {
			using enum mass_scheme_t;

			switch (p.channel.mass_scheme) {
			case zero_mass:
			case slow_rescaling:
			case slow_rescaling_modified:
				return sidis::coefficients::_impl::f3_hw::qg_nlo::lq_lq(p);
			case sacot_unsubtracted:
			case sacot_subtracted:
			case sacot_unsubtracted_modified:
			case sacot_subtracted_modified:
				return sidis::coefficients::_impl::f3_hw::qg_nlo::lq_hq(p);
			}
			std::unreachable();
		}

		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}
}

}

#endif