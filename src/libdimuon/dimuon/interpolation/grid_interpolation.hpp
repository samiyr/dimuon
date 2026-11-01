#ifndef DM_INTERPOLATION_GRID_INTERPOLATION_H
#define DM_INTERPOLATION_GRID_INTERPOLATION_H

#include <vector>
#include <algorithm>
#include <memory>
#include <optional>
#include <utility>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

namespace dm {

/// Wrapper for loading a grid file and interpolating it using GSL
struct grid_interpolation {
	/// Delegating constructor, see `grid_interpolation(std::vector<double>, std::vector<double>, const gsl_interp_type *)`
	explicit grid_interpolation(
		const std::pair<std::vector<double>, std::vector<double>> &grid_contents,
		const gsl_interp_type *interpolation_type
	) : grid_interpolation{std::get<0>(grid_contents), std::get<1>(grid_contents), interpolation_type} {}

	/// @brief Constructs the interpolator from the given data
	/// @param[in] grid_points List of grid points (i.e. the input points)
	/// @param[in] grid_values List of grid values (i.e. the output points, obtained from \p grid_points)
	/// @param[in] interpolation_type GSL interpolation type
	explicit grid_interpolation(
		std::vector<double> grid_points, 
		std::vector<double> grid_values, 
		const gsl_interp_type *interpolation_type
	) : spline{gsl_spline_alloc(interpolation_type, grid_points.size()), gsl_spline_free},
	accelerator{gsl_interp_accel_alloc(), gsl_interp_accel_free} {
		gsl_spline_init(spline.get(), grid_points.data(), grid_values.data(), grid_points.size());
	}

	/// @brief Returns the interpolated value at \p x
	/// @details Bounds checking is handled only by GSL. If an out-of-bounds input value is given, GSL will throw.
	double operator()(const double x) const noexcept {
		return gsl_spline_eval(spline.get(), x, accelerator.get());
	}

	/// Returns the interpolated value at \p x if it's in range, otherwise `std::nullopt`
	std::optional<double> try_interpolate(const double x) const noexcept {
		double value;
		if (gsl_spline_eval_e(spline.get(), x, accelerator.get(), &value) == GSL_EDOM) {
			return std::nullopt;
		}
		return value;
	}

	/// Returns the first and the last grid point
	std::pair<double, double> range() const noexcept {
		return {spline->interp->xmin, spline->interp->xmax};
	}

	/// Returns the interpolated value at \p x if it's in range, otherwise the value at the closest grid edge
	double interpolate_clamped(const double x) const noexcept {
		const double min = spline->interp->xmin;
		const double max = spline->interp->xmax;
		return operator()(std::clamp(x, min, max));
	}

private:
	std::shared_ptr<gsl_spline> spline;
	std::shared_ptr<gsl_interp_accel> accelerator;
};

}

#endif