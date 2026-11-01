#ifndef DM_COMMON_PARTONIC_CHANNEL_GROUP_H
#define DM_COMMON_PARTONIC_CHANNEL_GROUP_H

#include <cstddef>
#include <utility>
#include <vector>

namespace dm {

/// Position of a physical channel among the physical channels by which an evaluator (`dis::evaluator` or `sidis::evaluator`) indexes its partonic channels
using channel_index = std::size_t;

/// A group of partonic channels with the same coefficient functions (the same coefficient channel, see `dis::coefficient_channel` and `sidis::coefficient_channel`),
/// each with the position of its physical channel (see `channel_index`). Iteration gives the channels as `indexed_channel`s, in the order in which they were added.
/// @tparam partonic_channel_t Type of the partonic channels (`dis::partonic_channel` or `sidis::partonic_channel`)
template <typename partonic_channel_t>
class partonic_channel_group {
public:
	/// A partonic channel, with the position of its physical channel among the physical channels given to the evaluator constructor
	struct indexed_channel {
		partonic_channel_t channel;
		channel_index physical_channel_index;
	};

	/// Adds \p channel to the end of the group
	void push_back(indexed_channel channel) { channels.push_back(std::move(channel)); }

	/// The channels, in the order in which they were added
	decltype(auto) begin(this auto &&self) noexcept { return std::forward<decltype(self)>(self).channels.begin(); }
	/// The channels, in the order in which they were added
	decltype(auto) end(this auto &&self) noexcept { return std::forward<decltype(self)>(self).channels.end(); }

private:
	/// The channels, in the order in which they were added
	std::vector<indexed_channel> channels;
};

}

#endif
