#ifndef DM_DIS_PHYSICAL_CHANNEL_H
#define DM_DIS_PHYSICAL_CHANNEL_H

#include <compare> // IWYU pragma: keep (needed by the defaulted operator<=>)

#include "dimuon/common/flavor.hpp"
#include "dimuon/dis/partonic_channel.hpp"

namespace dm {

namespace dis {
	/// A physical channel, i.e. the flavor of the initial-state parton of a partonic channel.
	/// All partonic channels of a physical channel multiply the same PDF,
	/// so the physical channels carry all the information that a channel decomposition can give.
	struct physical_channel {
		/// Flavor of the initial-state parton, connected to the PDF
		flavor::type initial_flavor;

		/// Constructs the physical channel of \p channel
		explicit physical_channel(const partonic_channel &channel) noexcept : initial_flavor{channel.initial_line.flavor} {}

		/// Default the comparison operators to get member-wise ordering
		auto operator<=>(const physical_channel &) const = default;
	};
}

}

#endif
