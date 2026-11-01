#ifndef DM_COMMON_HEAVY_QUARK_THRESHOLDS_H
#define DM_COMMON_HEAVY_QUARK_THRESHOLDS_H

namespace dm {

/// Threshold of the heavy quarks, which determines the number of active flavors (relevant at NNLO).
/// These values are separate from the quark masses (see `parton`), although they default to the
/// quark masses of the PDF set.
struct heavy_quark_thresholds {
	/// Charm-quark threshold (GeV)
	double charm_threshold = 0.0;
	/// Bottom-quark threshold (GeV)
	double bottom_threshold = 0.0;
	/// Top-quark threshold (GeV)
	double top_threshold = 0.0;

	bool operator==(const heavy_quark_thresholds &) const = default;
};

}

#endif
