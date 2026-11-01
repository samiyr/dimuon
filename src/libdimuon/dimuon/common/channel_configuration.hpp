#ifndef DM_COMMON_CHANNEL_CONFIGURATION_H
#define DM_COMMON_CHANNEL_CONFIGURATION_H

#include "dimuon/common/mass_scheme.hpp"

namespace dm {

/// Configuration for an individual partonic channel
struct channel_configuration {
	/// Mass scheme used for the channel
	mass_scheme_t scheme;
};

}

#endif