import pandas as pd
import numpy as np

from common import *

### PDF and FF grid generation

def generate_slice(xfxQ2, flavor:int, knots:list[float], scales2:list[float]) -> float_array:
	return np.array([xfxQ2(flavor, knot, scale2) for scale2 in scales2 for knot in knots])

def pdf_slice_at_yE(pdf_xfxQ2, y:float, E:float, flavor:int, acc:accelerator) -> float_array:
	return generate_slice(pdf_xfxQ2, flavor, acc.x_knots, acc.factorization_scales[y, E])

def ff_slice_at_yE(ff_xfxQ2, y:float, E:float, flavor:int, acc:accelerator) -> float_array:
	return generate_slice(ff_xfxQ2, flavor, acc.z_knots, acc.fragmentation_scales[y, E])

### Convolutions

def convolve_over_flavors(df:pd.DataFrame, columns:list[str], pdf_xfxQ2, ff_xfxQ2, acc:accelerator) -> float_array:
	initial_flavors = acc.initial_flavors
	fragmentation_flavors = acc.fragmentation_flavors
	ys = acc.unique_ys
	Es = acc.unique_Es

	x_count = len(acc.unique_xs)
	y_count = len(ys)
	E_count = len(Es)
	xk_count = len(acc.x_knot_indices)
	zk_count = len(acc.z_knot_indices)
	if_count = len(initial_flavors)
	ff_count = len(fragmentation_flavors)

	pdf_matrix = np.array([pdf_slice_at_yE(pdf_xfxQ2, y, E, i_f, acc) for y in ys for E in Es for i_f in initial_flavors]).reshape(y_count, E_count, if_count, x_count, xk_count)
	ff_matrix = np.array([ff_slice_at_yE(ff_xfxQ2, y, E, f_f, acc) for y in ys for E in Es for f_f in fragmentation_flavors]).reshape(y_count, E_count, ff_count, x_count, zk_count)

	def matrix_product(column:str):
		data_matrix = slice_to_array(df, column).reshape(y_count, E_count, if_count, ff_count, x_count, xk_count, zk_count)
		return np.einsum("ijklmnp,ijkmn,ijlmp->ijm", data_matrix, pdf_matrix, ff_matrix, optimize=True)
	
	return np.array([matrix_product(column) for column in columns])

def convolve_over_resonances(grids:list[interpolation_grid], columns:list[str], pdf_xfxQ2, ff_xfxQ2s, ff_multipliers:list[float]) -> float_array:
	convolutions = np.array([convolve_over_flavors(grid.df, columns, pdf_xfxQ2, ff_xfxQ2, grid.acc) * mul for grid, ff_xfxQ2, mul in zip(grids, ff_xfxQ2s, ff_multipliers)])
	return np.nansum(convolutions, axis=0)

def convolve_df(
		grids:list[interpolation_grid],
		pdf_alphasQ2, pdf_xfxQ2, ff_xfxQ2s,
		ff_multipliers:list[float] = [1, 1, 1, 1.14]
	) -> pd.DataFrame:

	acc = grids[0].acc

	alpha_s = np.array([np.array([pdf_alphasQ2(renormalization_scale) for renormalization_scale in acc.renormalization_scales[y, E]]) for y in acc.unique_ys for E in acc.unique_Es]).reshape(len(acc.unique_ys), len(acc.unique_Es), -1)
	coeff = alpha_s / (2.0 * np.pi)

	convolutions = convolve_over_resonances(grids, ["LO", "NLO_correction", "NNLO_correction"], pdf_xfxQ2, ff_xfxQ2s, ff_multipliers)

	lo = convolutions[0].ravel()
	nlo_correction = coeff * convolutions[1]
	nnlo_correction = (coeff**2) * convolutions[2]

	nlo = lo + nlo_correction.ravel()
	nnlo = nlo + nnlo_correction.ravel()

	xs = acc.xs
	ys = acc.ys
	Es = acc.Es
	Q2s = acc.Q2s

	output_df = pd.DataFrame({
		"x": xs,
		"y": ys,
		"E": Es,
		"Q2": Q2s,
		"LO": lo,
		"NLO": nlo,
		"NNLO": nnlo
	})

	return output_df