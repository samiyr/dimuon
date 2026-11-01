#ifndef DM_SIDIS_COEFFICIENTS_F3_GC_H
#define DM_SIDIS_COEFFICIENTS_F3_GC_H

#include <utility>
#include <cmath> // IWYU pragma: keep (needed by the #include'd gc coefficient files)

#include "dimuon/sidis/coefficients/evaluation_parameters.hpp"
#include "dimuon/utility/math.hpp" // IWYU pragma: keep (needed by the #include'd gc coefficient files)

namespace dm {

/// The gc coefficient functions of F3. They are valid only for zeta >= z: they vanish below, where `sidis::evaluator::evaluate_observable` does not evaluate them.
namespace sidis::coefficients::f3_gc {
	namespace lo {
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/lo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C0q2qM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}

		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &) noexcept { return 0.0; }
		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}

	namespace nlo {
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nlo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C1q2qM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}

		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nlo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C1q2gM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}

		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nlo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C1g2qM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}

		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}
	
	namespace nnlo {
		inline vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nnlo_decl.hpp"

				switch(p.channel.c_symmetric_type) {
					case partonic_channel::channel_type::beta_alpha:
					case partonic_channel::channel_type::alphabar_betabar:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qM.hpp"
					}
					case partonic_channel::channel_type::qq_Fcon:
					case partonic_channel::channel_type::qbar_qbar_Fcon:
						return 0.0;
					case partonic_channel::channel_type::alpha_alpha_Fcon:
					case partonic_channel::channel_type::betabar_betabar_Fcon:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qMFcon1.hpp"
					}
					case partonic_channel::channel_type::betabar_alpha:
					case partonic_channel::channel_type::alpha_betabar:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qbM.hpp"
					}
					case partonic_channel::channel_type::alphabar_alpha_Fcon:
					case partonic_channel::channel_type::beta_betabar_Fcon:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qbMFcon.hpp"
					}
					case partonic_channel::channel_type::qprime_alpha:
					case partonic_channel::channel_type::qbarprime_alpha:
					case partonic_channel::channel_type::qprime_betabar:
					case partonic_channel::channel_type::qbarprime_betabar:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qpM1.hpp"
					}
					case partonic_channel::channel_type::betaprime_q:
					case partonic_channel::channel_type::alphabarprime_q:
					case partonic_channel::channel_type::betaprime_qbar:
					case partonic_channel::channel_type::alphabarprime_qbar:
					{
						#include "dimuon/sidis/coefficients/gc/F3C2q2qpM2.hpp"
					}
					case partonic_channel::channel_type::g_alpha:
					case partonic_channel::channel_type::g_betabar:
					case partonic_channel::channel_type::beta_g:
					case partonic_channel::channel_type::alphabar_g:
					case partonic_channel::channel_type::gg:
					{
						std::unreachable();
					}
				}
				std::unreachable();
			}() / p.z / p.pdf_momentum_fraction;
		}
		inline vectorization::distribution_product quark_to_gluon(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nnlo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C2q2gM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}
		inline vectorization::distribution_product gluon_to_quark(const sidis::evaluation_parameters &p) noexcept {
			return [&p]() -> vectorization::distribution_product {
				#include "dimuon/sidis/coefficients/gc/nnlo_decl.hpp"
				#include "dimuon/sidis/coefficients/gc/F3C2g2qM.hpp"
			}() / p.z / p.pdf_momentum_fraction;
		}
		inline vectorization::distribution_product gluon_to_gluon(const sidis::evaluation_parameters &) noexcept { return 0.0; }
	}
}

}

#endif