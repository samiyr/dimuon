#ifndef DM_COMMON_IS_OBSERVABLE_H
#define DM_COMMON_IS_OBSERVABLE_H

#include "dimuon/common/trf_kinematics.hpp"
#include "dimuon/common/scattering_process.hpp"
#include "dimuon/common/cross_section_observable.hpp"

namespace dm {

/// An observable-type functor should be a callable type that takes the kinematical point, scattering process,
/// and the values of F1, F2, and F3, and returns a number.
/// Mathematically, an observable has to be linear in F1, F2, and F3, with no constant term.
template <typename T>
concept is_observable = requires(T observable, const trf_kinematics &kinematics, const scattering_process &process, const double f1, const double f2, const double f3) {
	{ observable(kinematics, process, f1, f2, f3) } -> std::convertible_to<double>;
};

static_assert(is_observable<cross_section_observable>, "Type `cross_section_observable` doesn't satisfy concept `is_observable`.");

}

#endif