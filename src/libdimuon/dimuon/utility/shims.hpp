#ifndef DM_UTILITY_SHIMS_H
#define DM_UTILITY_SHIMS_H

#include <version> // IWYU pragma: keep (defines the __cpp_lib_* feature-test macros)

#ifdef __cpp_lib_inplace_vector

#include <inplace_vector>

namespace dm::shims {
	using std::inplace_vector;
}

#else

#include <inplace_vector/inplace_vector.hpp>

namespace dm::shims {
	using beman::inplace_vector::inplace_vector;
}

#endif

#ifdef __cpp_lib_mdspan

#include <mdspan>

namespace dm::shims {
	using std::mdspan;
	using std::dextents;
	using std::layout_right;
}

#else

#include <mdspan/mdspan/mdspan.hpp> // IWYU pragma: keep

namespace dm::shims {
	using Kokkos::mdspan;
	using Kokkos::dextents;
	using Kokkos::layout_right;
}

#endif

#ifdef __cpp_lib_submdspan

#include <mdspan>

namespace dm::shims {
	using std::submdspan;
	using std::full_extent;
}

#else

#include <mdspan/mdspan/mdspan.hpp> // IWYU pragma: keep

namespace dm::shims {
	using Kokkos::submdspan;
	using Kokkos::full_extent;
}

#endif

// `slice_view` implementation is unconditional because there's no feature-test macro for `std::views::slice` yet
#include <slice_view/slice_view.hpp> // IWYU pragma: keep

namespace dm::shims::views {
	using namespace beman::slice_view::views;
}

#endif