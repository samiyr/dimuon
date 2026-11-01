#ifndef DM_PDF_INTERFACE_VECTOR_KERNEL_INTERFACE_H
#define DM_PDF_INTERFACE_VECTOR_KERNEL_INTERFACE_H

#include <cstddef>
#include <utility>

#include "dimuon/common/flavor.hpp"

#include "dimuon/interpolation/linear_interpolation_basis.hpp"

#include "dimuon/utility/sparse_vector.hpp"

namespace dm {

/// Interface for computing linear interpolation kernel values in a vectorized form
class vector_kernel_interface {
public:
	/// The type of output values of all kernels at a given point
	using values = sparse_vector<double, linear_interpolation_basis::maximum_kernels_per_point>;

	/// @brief Constructs the interface from the kernels
	/// @param[in] kernels_ The interpolation kernels
	explicit vector_kernel_interface(linear_interpolation_basis kernels_) noexcept : kernels{std::move(kernels_)} {}

	/// Dimension of the output vector
	std::size_t dimension() const noexcept { return kernels.size(); }

	/// Returns the values of the kernels at \p x, which are the same for all flavors and scales
	values xfxQ2(const flavor::type, const double x, const double Q2) const { return xfxQ2(x, Q2); }

	/// Returns the values of the kernels at \p x, which are the same for all flavors and scales
	values xfxQ2(const double x, const double) const {
		values kernel_values;
		for (const linear_interpolation_basis::kernel &kernel : kernels.nonzero_at(x)) {
			kernel_values[kernel.index] = kernel(x);
		}
		return kernel_values;
	}

	/// Kernels have no `alpha_s`, so this function always returns zero
	double alpha_s(const double) const noexcept { return 0.0; }

private:
	/// The interpolation kernels
	linear_interpolation_basis kernels;
};

}

#endif
