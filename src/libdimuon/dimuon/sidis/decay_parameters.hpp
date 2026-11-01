#ifndef DM_SIDIS_DECAY_PARAMETERS_H
#define DM_SIDIS_DECAY_PARAMETERS_H

#include <vector>

#include "dimuon/decay/decay_instance.hpp"
#include "dimuon/decay/decay_function.hpp"

namespace dm {

/// @brief Struct for containing parameters related to the decay process
/// @tparam decay_function_t Type of decay function object, used to get the value of the decay function. Must satisfy the `is_decay_function` concept.
template <is_decay_function decay_function_t>
struct decay_parameters {
	/// Array of decay instances, one for each decay channel
	std::vector<decay_instance<decay_function_t>> decay_instances;
};

}

#endif