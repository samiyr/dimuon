#ifndef DM_EXE_RUN_MODE_H
#define DM_EXE_RUN_MODE_H

/// Run mode for the observable calculation
enum class run_mode {
	/// Uses a PDF and FFs and computes the cross section directly
	direct,
	/// Uses PDF and FF interpolation kernels
	kernel
};

#endif