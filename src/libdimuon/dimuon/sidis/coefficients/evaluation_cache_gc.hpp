#ifndef DM_SIDIS_COEFFICIENTS_EVALUATION_CACHE_GC_H
#define DM_SIDIS_COEFFICIENTS_EVALUATION_CACHE_GC_H

#include <cmath>
#include <numbers> // IWYU pragma: keep (needed by the #include'd gc cache files)

#include "dimuon/common/heavy_quark_thresholds.hpp"
#include "dimuon/common/perturbative_order.hpp"
#include "dimuon/common/constants.hpp"
#include "dimuon/utility/math.hpp"

namespace dm {

namespace sidis {
	/// Values of the gc coefficient functions that do not depend on the channel, only on the integration point and the kinematical point.
	/// They are computed once per integration point at construction, instead of in each coefficient function, which reads them as `p.c.<order>.<name>`
	/// (see `gc/{lo,nlo,nnlo}_decl.hpp`). The values of each order are generated with the coefficients (see `gc/{lo,nlo,nnlo}_cache.hpp`).
	class evaluation_cache_gc {
	public:
		/// The quantities that the values depend on. `xi`, `zeta`, and `z` have the names of the fields of `sidis::evaluation_parameters`, since the generated expressions read them from both.
		struct inputs {
			/// The PDF-side convolution variable
			double xi;
			/// The FF-side convolution variable
			double zeta;
			/// The FF-side momentum fraction
			double z;
			/// The renormalization scale at the current kinematical point
			double renormalization_scale;
			/// The factorization scale at the current kinematical point
			double factorization_scale;
			/// The fragmentation scale at the current kinematical point
			double fragmentation_scale;
			/// The virtuality of the current kinematical point
			double Q2;
			/// The thresholds of the heavy quarks (GeV), which determine the number of active flavors
			heavy_quark_thresholds thresholds;
		};

		/// The LO values are constants, so they are set at construction
		struct lo_cache {
			lo_cache() noexcept {
				#include "dimuon/sidis/coefficients/gc/lo_cache.hpp"
			}

			#include "dimuon/sidis/coefficients/gc/lo_vars.hpp"
		};

		struct nlo_cache {
			void compute(const sidis::evaluation_cache_gc::inputs &p, const evaluation_cache_gc *c) noexcept {
				using constants::NC;

				const double LMUF = c->LMUF;
				const double LMUA = c->LMUA;

				#include "dimuon/sidis/coefficients/gc/nlo_cache.hpp"
			}

			#include "dimuon/sidis/coefficients/gc/nlo_vars.hpp"
		};

		struct nnlo_cache {
			void compute(const sidis::evaluation_cache_gc::inputs &p, const evaluation_cache_gc *c) noexcept {
				using constants::NC;

				const double LMUR = c->LMUR;
				const double LMUF = c->LMUF;
				const double LMUA = c->LMUA;

				const double NF = c->NF;
				const double poly2 = c->poly2;
				const double sqrt1 = c->sqrt1;
				const double sqrt2 = c->sqrt2;
				const double sqrt3 = c->sqrt3;

				enum class sidis_nnlo_region_t { r1, r2, t1, t2, u1, u2, u3, u4 };

				constexpr sidis_nnlo_region_t r1 = sidis_nnlo_region_t::r1;
				constexpr sidis_nnlo_region_t r2 = sidis_nnlo_region_t::r2;

				constexpr sidis_nnlo_region_t t1 = sidis_nnlo_region_t::t1;
				constexpr sidis_nnlo_region_t t2 = sidis_nnlo_region_t::t2;

				constexpr sidis_nnlo_region_t u1 = sidis_nnlo_region_t::u1;
				constexpr sidis_nnlo_region_t u2 = sidis_nnlo_region_t::u2;
				constexpr sidis_nnlo_region_t u3 = sidis_nnlo_region_t::u3;
				constexpr sidis_nnlo_region_t u4 = sidis_nnlo_region_t::u4;

				// The theta function of each region of the (xi, zeta) plane in which the NNLO coefficient functions have a different form, as defined
				// in arXiv:2201.06982, eqs. (5.9), (5.12), and (5.16), with x and z there being xi and zeta here: R1 and R2 are the two sides of
				// zeta = 1 - xi, T1 and T2 the two sides of zeta = xi, and U1-U4 their intersections
				const auto T = [&p](const sidis_nnlo_region_t region) -> double {
					switch(region) {
						case sidis_nnlo_region_t::r1: return p.zeta < 1.0 - p.xi;
						case sidis_nnlo_region_t::r2: return p.zeta > 1.0 - p.xi;

						case sidis_nnlo_region_t::t1: return p.zeta > p.xi;
						case sidis_nnlo_region_t::t2: return p.zeta < p.xi;

						case sidis_nnlo_region_t::u1: return p.zeta < 1.0 - p.xi && p.zeta < p.xi;
						case sidis_nnlo_region_t::u2: return p.zeta > 1.0 - p.xi && p.zeta < p.xi;
						case sidis_nnlo_region_t::u3: return p.zeta < 1.0 - p.xi && p.zeta > p.xi;
						case sidis_nnlo_region_t::u4: return p.zeta > 1.0 - p.xi && p.zeta > p.xi;
					}
					std::unreachable();
				};

				#include "dimuon/sidis/coefficients/gc/nnlo_cache.hpp"
			}

			#include "dimuon/sidis/coefficients/gc/nnlo_vars.hpp"
		};

		/// An empty cache, for coefficient functions that don't read it (the hw implementation). Its values are uninitialized, apart from the constant LO values.
		evaluation_cache_gc() = default;

		/// Computes the values necessary for the evaluation of the coefficient functions up to order \p pto for inputs \p p. The other values are uninitialized.
		/// - LO: the values are constants, which `lo_cache` sets at its construction.
		/// - NLO and above: this function computes the shared transcendental functions 
		/// (logs at NLO, also polylogs, arctangents, and inverse tangent integrals at NNLO) declared in `sidis/coefficients/gc/{pto}_vars.hpp`, and the
		/// three scales (renormalization, factorization, and fragmentation).
		/// - NNLO only: this function computes certain other variables shared among the NNLO coefficients.
		explicit evaluation_cache_gc(const inputs &p, const perturbative_order pto) noexcept {
			if (pto >= perturbative_order::nlo) {
				// Renormalization scale
				LMUR = std::log(p.renormalization_scale / p.Q2);
				// Factorization scale
				LMUF = std::log(p.factorization_scale / p.Q2);
				// Fragmentation scale
				LMUA = std::log(p.fragmentation_scale / p.Q2);

				nlo.compute(p, this);
			}

			if (pto >= perturbative_order::nnlo) {
				poly2 = 1.0 + 2.0 * p.xi + math::power<2>(p.xi) - 4.0 * p.xi * p.zeta;
				sqrt1 = std::sqrt(1.0 - 2.0 * p.zeta + math::power<2>(p.zeta) + 4.0 * p.xi * p.zeta);
				sqrt2 = std::sqrt(poly2);
				sqrt3 = std::sqrt(p.xi / p.zeta);

				// Compute the number of active flavors
				NF = [&p]() -> double {
					const double scale = p.factorization_scale;

					double nf = 3.0;
					if (scale >= math::power<2>(p.thresholds.charm_threshold)) { nf += 1.0; }
					if (scale >= math::power<2>(p.thresholds.bottom_threshold)) { nf += 1.0; }
					if (scale >= math::power<2>(p.thresholds.top_threshold)) { nf += 1.0; }

					return nf;
				}();

				nnlo.compute(p, this);
			}
		}

		lo_cache lo;
		nlo_cache nlo;
		nnlo_cache nnlo;

		double LMUR;
		double LMUF;
		double LMUA;

		double NF;
		double poly2;
		double sqrt1;
		double sqrt2;
		double sqrt3;
	};
}

}

#endif
