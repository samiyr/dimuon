#ifndef DM_EXE_SCALE_DEPENDENCE_MODE_H
#define DM_EXE_SCALE_DEPENDENCE_MODE_H

/// Scale dependence mode: the renormalization, factorization, and fragmentation scales (before applying multipliers and lower bounds)
enum class scale_dependence_mode {
	/// mu^2 = Q^2
	Q2,
	/// mu^2 = Q^2 + mc^2
	Q2mc2
};

#endif
