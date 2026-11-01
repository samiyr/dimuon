#ifndef DM_COUPLING_IS_COUPLING_H
#define DM_COUPLING_IS_COUPLING_H

#include <concepts>

#include "dimuon/common/flavor.hpp"

#include "dimuon/coupling/ckm_matrix.hpp"
#include "dimuon/coupling/electromagnetic_coupling.hpp"

namespace dm {

/// A coupling-type functor should be a callable type that takes two flavors and returns a number
template <typename T>
concept is_coupling = requires(T coupling, const flavor::type flavor1, const flavor::type flavor2) {
	{ coupling(flavor1, flavor2) } -> std::convertible_to<double>;
};

static_assert(is_coupling<ckm_matrix>, "Type `ckm_matrix` doesn't satisfy concept `is_coupling`");
static_assert(is_coupling<electromagnetic_coupling>, "Type `electromagnetic_coupling` doesn't satisfy concept `is_coupling`");

}

#endif