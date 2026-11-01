#ifndef DM_PDF_EXTRAPOLATOR_ZERO_EXTRAPOLATOR_H
#define DM_PDF_EXTRAPOLATOR_ZERO_EXTRAPOLATOR_H

#include "LHAPDF/Extrapolator.h"

namespace dm {

/// LHAPDF extrapolator that's zero outside the grid
class zero_extrapolator : public LHAPDF::Extrapolator {
public:
	double extrapolateXQ2(int, double, double) const noexcept override {
		return 0.0;
	}
};

}

#endif