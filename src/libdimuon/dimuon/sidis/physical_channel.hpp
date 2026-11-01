#ifndef DM_SIDIS_PHYSICAL_CHANNEL_H
#define DM_SIDIS_PHYSICAL_CHANNEL_H

#include <compare> // IWYU pragma: keep (needed by the defaulted operator<=>)

#include "dimuon/common/flavor.hpp"
#include "dimuon/sidis/partonic_channel.hpp"

namespace dm {

namespace sidis {
	/// A physical channel, i.e. the flavors of the initial-state and fragmenting partons of a partonic channel.
	/// All partonic channels of a physical channel multiply the same product of a PDF and an FF,
	/// so the physical channels carry all the information that a channel decomposition can give.
	struct physical_channel {
		/// Flavor of the initial-state parton, connected to the PDF
		flavor::type initial_flavor;
		/// Flavor of the fragmenting parton, connected to the FF
		flavor::type fragmentation_flavor;

		/// Constructs the physical channel of \p channel
		explicit physical_channel(const partonic_channel &channel) noexcept :
			initial_flavor{channel.initial_line.flavor}, fragmentation_flavor{channel.fragmentation_line.flavor} {}

		/// Default the comparison operators to get member-wise ordering, first by the initial and then by the fragmenting flavor
		auto operator<=>(const physical_channel &) const = default;
	};
}

}

#endif
