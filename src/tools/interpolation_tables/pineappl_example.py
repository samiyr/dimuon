import pineappl
import numpy as np
import pandas as pd
import lhapdf
from pathlib import Path

def pineappl_cross_section(pdf, ffs, grids):
	space_like_type = pineappl.convolutions.ConvType(polarized=False, time_like=False)
	time_like_type = pineappl.convolutions.ConvType(polarized=False, time_like=True)

	types = [
		pineappl.convolutions.Conv(convolution_types=space_like_type, pid=2212),
		pineappl.convolutions.Conv(convolution_types=time_like_type, pid=2212)
	]

	def get_kinematic_rows(grid):
		return np.array([[r[0] for r in row] for row in grid.bin_limits()])

	def convolve(grids, ffs, order_mask):
		convolutions = [grid.convolve(
			pdg_convs=types,
			xfxs=[pdf.xfxQ2, ff],
			alphas=pdf.alphasQ2,
			order_mask=order_mask
		) for ff, grid in zip(ffs, grids)]
		
		return np.sum(convolutions, axis=0)

	all_kinematic_rows = [get_kinematic_rows(grid) for grid in grids]
	kinematic_rows = pd.DataFrame(all_kinematic_rows[0], columns=["x", "y", "E"])

	if not all((x == all_kinematic_rows[0]).all() for x in all_kinematic_rows):
		raise Exception("The grids must have identical kinematics")

	lo = convolve(grids, ffs, [True, False, False])
	nlo = convolve(grids, ffs, [True, True, False])
	nnlo = convolve(grids, ffs, [True, True, True])

	pto_data = pd.DataFrame({"LO" : lo, "NLO" : nlo, "NNLO" : nnlo})

	return pd.concat([kinematic_rows, pto_data], axis=1)

if __name__ == '__main__':
	pdf = lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Fe56")

	ffs = [
		lhapdf.mkPDF(name)
		for name in [
			"kkks08_opal_d0___mas",
			"kkks08_opal_d+___mas",
			"bkk05_D3_d_s_nlo",
			"bkk05_D3_lambda_c_nlo"
		]
	]

	# Halve the light-quark flavors when using KKKS08 / BKK05 FFs
	multipliers = {-6 : 1.0, -5: 1.0, -4: 1.0, -3: 0.5, -2: 0.5, -1: 0.5, 0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5, 4: 1.0, 5: 1.0, 6: 1.0}

	def bind_ff(ff, multipliers):
		return lambda f, x, Q2: ff.xfxQ2(f, x, Q2) * multipliers[f]

	ff_funcs = [bind_ff(ff, multipliers) for ff in ffs]

	grids = [pineappl.grid.Grid.read(Path(f"path/to/eppslike/grid/{res}.pineappl.lz4")) for res in ["D0", "Dp", "Ds", "LambdaC"]]

	print(pineappl_cross_section(pdf, ff_funcs, grids))