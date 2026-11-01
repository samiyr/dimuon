#ifndef DM_PDF_EXTRAPOLATOR_FREEZE_EXTRAPOLATOR_H
#define DM_PDF_EXTRAPOLATOR_FREEZE_EXTRAPOLATOR_H

#include <algorithm>

#include "LHAPDF/Extrapolator.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/KnotArray.h"

namespace dm {

/// LHAPDF extrapolator that freezes at the edges of the grid: `x` and `Q^2` outside the grid are clamped to the first or last knot
class freeze_extrapolator : public LHAPDF::Extrapolator {
public:

	double extrapolateXQ2(int id, double x, double q2) const override {
		// Clamped to the outermost knots exactly: LHAPDF's interpolator doesn't check the range, and reads outside its arrays below the first knot
		const LHAPDF::KnotArray &knots = pdf().knotarray();
		const double clamped_x = std::clamp(x, knots.xs().front(), knots.xs().back());
		const double clamped_Q2 = std::clamp(q2, knots.q2s().front(), knots.q2s().back());
		return pdf().interpolator().interpolateXQ2(id, clamped_x, clamped_Q2);
	}
};

}

#endif