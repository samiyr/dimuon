#ifndef DM_EXE_INTERPOLATION_KERNEL_H
#define DM_EXE_INTERPOLATION_KERNEL_H

/// Interpolation kernel types
enum class interpolation_kernel {
	/// Piecewise-linear kernels (see `linear_interpolation_basis`)
	linear,
	/// Global Lagrange kernels (see `lagrange_interpolation_basis`)
	global_lagrange,
	/// Composite Lagrange kernels (i.e. `composite_interpolation_basis<lagrange_interpolation_basis>`),
	/// where Lagrange kernels are confined to blocks
	composite_lagrange
};

#endif
