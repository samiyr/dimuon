#ifndef DM_COMMON_STRUCTURE_FUNCTION_H
#define DM_COMMON_STRUCTURE_FUNCTION_H

#include <utility>

#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/common/scattering_process.hpp"

namespace dm {

/// Structure function type (F1, F2, F3)
enum class structure_function_t { f1, f2, f3 };

namespace structure_function {
	struct basis_123_representation {
		double f1;
		double f2;
		double f3;
	};
	struct basis_tl3_representation {
		double ft;
		double fl;
		double f3;

		constexpr explicit basis_tl3_representation(const basis_123_representation &b123, double pdf_momentum_fraction) noexcept
		: ft(2.0 * b123.f1), fl(b123.f2 / pdf_momentum_fraction - 2.0 * b123.f1), f3(b123.f3) {}
	};

	template <structure_function_t sf>
	struct observable {
		constexpr double operator()(const trf_kinematics &, const scattering_process &, const double f1, const double f2, const double f3) const noexcept {
			if constexpr (sf == structure_function_t::f1) { return f1; }
			else if constexpr (sf == structure_function_t::f2) { return f2; }
			else if constexpr (sf == structure_function_t::f3) { return f3; }
			std::unreachable();
		}
	};
}

}

#endif