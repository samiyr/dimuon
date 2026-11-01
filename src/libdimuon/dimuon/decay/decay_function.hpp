#ifndef DM_DECAY_FUNCTION_H
#define DM_DECAY_FUNCTION_H

#include <concepts>

#include "dimuon/common/particle.hpp"

namespace dm {

/// A decay function should take in the `x`, `z`, and `Q2` kinematical values and a `target` particle, and return a number
template <typename T>
concept is_decay_function = requires(
	T decay_function, const double x, const double z, const double Q2, const particle &target) {
	{ decay_function(x, z, Q2, target) } -> std::convertible_to<double>;
};

}

#endif