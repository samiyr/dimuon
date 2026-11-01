import itertools
from pathlib import Path

import numpy as np

from common import interpolation_grid, import_gridfile

import pineappl

### PineAPPL conversion

def convert_grid_to_pineappl(input_grid:interpolation_grid, normalization:float=1.0) -> pineappl.grid.Grid:
	"""
	Converts the given input grid into a PineAPPL grid
	
	Parameters
	----------
	input_grid : interpolation_grid 
		The input grid object for a single resonance, obtained using e.g. `common.import_gridfile` or `common.import_grid`
	normalization : float [default: 1.0]
		Normalization constant multiplying all appropriate values in the output grid	

	Returns
	-------
	pineappl.grid.Grid
		PineAPPL grid object
	"""

	df = input_grid.df
	acc = input_grid.acc

	# Partonic x knots (corresponding to the interpolation of the PDF)
	x_knots = acc.x_knots
	# Partonic z knots (corresponding to the interpolation of the FF)
	z_knots = acc.z_knots

	# Set up perturbative orders and couplings
	pto_orders = list(range(0, 3))
	order_columns = ["LO", "NLO_correction", "NNLO_correction"]
	orders = [pineappl.boc.Order(qcd_pto, 0, 0, 0, 0) for qcd_pto in pto_orders]

	# Set up channels
	# Only the flavor combinations that are present in the grid, in the order of the full product of flavors
	present_pids = set(df[["if", "ff"]].drop_duplicates().itertuples(index=False, name=None))
	channel_pids = [
		(int(i_f), int(f_f))
		for i_f, f_f in itertools.product(acc.initial_flavors, acc.fragmentation_flavors)
		if (i_f, f_f) in present_pids
	]
	channels = [pineappl.boc.Channel([([i_f, f_f], 1.0)]) for i_f, f_f in channel_pids]

	# PDF convolution
	space_like_type = pineappl.convolutions.ConvType(polarized=False, time_like=False)
	# FF convolution
	time_like_type = pineappl.convolutions.ConvType(polarized=False, time_like=True)

	# Define convolveable objects (PDF and FF)
	convolutions = [
		pineappl.convolutions.Conv(convolution_types=space_like_type, pid=2212),
		pineappl.convolutions.Conv(convolution_types=time_like_type, pid=2212)
	]

	# Set up actual kinematics of the dimuon process
	xs = acc.unique_xs
	ys = acc.unique_ys
	Es = acc.unique_Es

	kinematical_points = list(itertools.product(Es, ys, xs))

	# Skeleton kinematics bins, to be filled properly later
	bins = len(kinematical_points)
	bin_limits = pineappl.boc.BinsWithFillLimits.from_fill_limits(
		fill_limits=list(map(float, range(0, bins + 1)))
	)

	# Set up "kinematical" variables
	# The actual kinematical (in terms of the dimuon process) variables
	# are defined above with the bins
	kinematics = [
		pineappl.boc.Kinematics.Scale(0), 	# renormalization scale
		pineappl.boc.Kinematics.Scale(1), 	# factorization scale
		pineappl.boc.Kinematics.Scale(2), 	# fragmentation scale
		pineappl.boc.Kinematics.X(0), 		# partonic x grid
		pineappl.boc.Kinematics.X(1)		# partonic z grid
	]
	
	# Set up scale functions
	scale_functions = pineappl.boc.Scales(
		# renormalization scale function, connected to the renormalization scale via the index 0
		ren=pineappl.boc.ScaleFuncForm.Scale(0),
		# factorization scale function, connected to the factorization scale via the index 1
		fac=pineappl.boc.ScaleFuncForm.Scale(1),	
		# fragmentation scale function, connected to the fragmentation scale via the index 2
		frg=pineappl.boc.ScaleFuncForm.Scale(2)
	)

	# Array of renormalization/factorization/fragmentation scales, obtained directly from the grid file
	renormalization_scales = np.array(list(acc.renormalization_scales.values()))
	factorization_scales = np.array(list(acc.factorization_scales.values()))
	fragmentation_scales = np.array(list(acc.fragmentation_scales.values()))

	# Configuration for the interpolation of scale values
	scale_nodes = 50
	scale_order = 3
	scale_reweight_meth = pineappl.interpolation.ReweightingMethod.NoReweight
	scale_map = pineappl.interpolation.MappingMethod.ApplGridH0
	scale_interpolation_meth = pineappl.interpolation.InterpolationMethod.Lagrange

	# Set up interpolation for the "kinematical" variables
	interpolations = [
		pineappl.interpolation.Interp(
			min=renormalization_scales.min(),
			max=renormalization_scales.max(),
			nodes=scale_nodes,
			order=scale_order,
			reweight_meth=scale_reweight_meth,
			map=scale_map,
			interpolation_meth=scale_interpolation_meth
		), # renormalization scale

		pineappl.interpolation.Interp(
			min=factorization_scales.min(),
			max=factorization_scales.max(),
			nodes=scale_nodes,
			order=scale_order,
			reweight_meth=scale_reweight_meth,
			map=scale_map,
			interpolation_meth=scale_interpolation_meth
		), # factorization scale

		pineappl.interpolation.Interp(
			min=fragmentation_scales.min(),
			max=fragmentation_scales.max(),
			nodes=scale_nodes,
			order=scale_order,
			reweight_meth=scale_reweight_meth,
			map=scale_map,
			interpolation_meth=scale_interpolation_meth
		), # fragmentation scale

		pineappl.interpolation.Interp(
			min=x_knots.min(),
			max=x_knots.max(),
			nodes=len(x_knots),
			order=4,
			reweight_meth=pineappl.interpolation.ReweightingMethod.ApplGridX,
			map=pineappl.interpolation.MappingMethod.ApplGridF2,
			interpolation_meth=pineappl.interpolation.InterpolationMethod.Lagrange
		), # partonic x knots

		pineappl.interpolation.Interp(
			min=z_knots.min(),
			max=z_knots.max(),
			nodes=len(z_knots),
			order=4,
			reweight_meth=pineappl.interpolation.ReweightingMethod.ApplGridX,
			map=pineappl.interpolation.MappingMethod.ApplGridF2,
			interpolation_meth=pineappl.interpolation.InterpolationMethod.Lagrange
		) # partonic z knots
	]

	# Set up PineAPPL grid
	# This is somewhat of a skeleton at this stage,
	# the actual data is filled next
	grid = pineappl.grid.Grid(
		pid_basis=pineappl.pids.PidBasis.Pdg,
		channels=channels,
		orders=orders,
		bins=bin_limits,
		convolutions=convolutions,
		interpolations=interpolations,
		kinematics=kinematics,
		scale_funcs=scale_functions
	)

	# Accumulation of actual kinematical variable bin limits
	limits = []

	# Enumerate over actual kinematical points
	for bin_idx, (E, y, x) in enumerate(kinematical_points):
		# Add to the kinematical bin limits
		# The limits are trivial, i.e. bin width is zero
		limits.append([(x, x), (y, y), (E, E)])
		
		# Slice of the full grid to the current kinematical point
		df_slice = df[np.logical_and.reduce([
				df["x"] == x,
				df["y"] == y,
				df["E"] == E
			]
		)]

		# Get the scale values of the current kinematical point from the grid
		# This just gets the first element, since the scales are functions of Q^2
		# which is fully determined from the triplet (x, y, E) - that is to say,
		# the scale values in df_slice are all the same and it doesn't matter which
		# one is taken
		renormalization_scale = df_slice["ren"].iloc[0]
		factorization_scale = df_slice["fact"].iloc[0]
		fragmentation_scale = df_slice["frag"].iloc[0]

		# Enumerate over the partonic channels (i_f = initial flavor, f_f = fragmentation flavor)
		for channel_idx, (i_f, f_f) in enumerate(channel_pids):
			# Another slice, now to the current partonic channel
			channel_slice = df_slice[np.logical_and(df_slice["if"] == i_f, df_slice["ff"] == f_f)]

			# A channel missing at this kinematical point (possible when merging grid files) is left as an empty subgrid
			if channel_slice.empty:
				continue

			# Enumerate over the perturbative orders
			# order_idx == 0: LO
			# order_idx == 1: NLO correction
			# order_idx == 2: NNLO correction
			for order_idx, order_column in enumerate(order_columns):
				# PineAPPL calculates higher-order corrections by multiplying by alpha_s^n
				# The grid that is to be converted expects multiplication by (alpha_s / 2π)^n,
				# hence this additional multiplicative factor 
				coeff = (1.0 / (2.0 * np.pi)) ** order_idx

				# At this point, the full grid can be reduced to a 2D matrix, populated with cross section
				# values corresponding to a pair of x knot and z knot values. The shape of this matrix is
				# (number of x knots) by (number of z knots)
				matrix_shape = len(x_knots), len(z_knots)

				# A matrix of x and z knot values, with the same shape as the actual data matrix. This is
				# just because of what PineAPPL expects f(x) and not x*f(x) for PDFs (and same for FFs)
				knot_array = (channel_slice["xk"] * channel_slice["zk"]).to_numpy().reshape(matrix_shape)

				# Data matrix
				matrix = normalization * coeff * channel_slice[order_column].to_numpy().reshape(matrix_shape) * knot_array

				# Add 3 trivial dimensions to the 2D data matrix to accomodate the addition of the scale values
				# New shape is (1, 1, 1, matrix_shape[0], matrix_shape[1])
				array = matrix[np.newaxis, np.newaxis, np.newaxis, ...]

				# 5D matrix with the scale values and partonic knots
				node_values = [np.array([renormalization_scale]), np.array([factorization_scale]), np.array([fragmentation_scale]), x_knots, z_knots]

				# Add subgrid to the PineAPPL grid
				subgrid = pineappl.subgrid.ImportSubgridV1(
					array=array,
					node_values=node_values
				)
				grid.set_subgrid(order_idx, bin_idx, channel_idx, subgrid.into())

	# Set up actual bin limits, replacing previous skeleton bins
	bin_configs = pineappl.boc.BinsWithFillLimits.from_limits_and_normalizations(
		limits=limits,
		normalizations=[1.0] * bins
	)
	grid.set_bwfl(bin_configs)

	grid.set_metadata("x1_label", "x")
	grid.set_metadata("x2_label", "y")
	grid.set_metadata("x3_label", "E")
	grid.set_metadata("x3_unit", "GeV")
	grid.set_metadata("y_label", "xs")

	grid.optimize()
	return grid

def convert_gridfile_to_pineappl(input_path:Path, output_directory:Path, resonances:dict[str, float]={"D0": 1.0, "Dp": 1.0, "Ds": 1.0, "LambdaC": 1.14}) -> None:
	"""
	Converts a single grid file to PineAPPL grids, one for each resonance, and writes them to file

	Parameters
	----------
	input_path:Path
		Path pointing to the grid file
	output_directory:Path
		Path pointing to the directory in which the PineAPPL grids will be written
	resonances:dict[str, float] [default: {"D0": 1.0, "Dp": 1.0, "Ds": 1.0, "LambdaC": 1.14}]
		Dictionary defining the resonance particle names (which must match the ones in the input grids) and their corresponding normalization constants,
		which will multiply all appropriate values in the output grids

	Notes
	-----
	The filenames for the output grids are of the format '<resonance>.pineappl.lz4', where <resonance> is the name of the resonance particle defined in `resonances`
	"""

	print(f"Importing {input_path}")
	# Missing flavor combinations are not needed, since they would only be removed as empty channels in `grid.optimize()`
	grids = {grid.acc.unique_resonances[0]: grid for grid in import_gridfile(input_path, fill_missing_channels=False)}

	# Resonances requested in `resonances` but missing from the grid file are skipped
	missing_resonances = [res for res in resonances.keys() if res not in grids]
	if missing_resonances:
		print(f"Warning: no data found for resonances {', '.join(missing_resonances)}, skipping them")

	output_directory.mkdir(parents=True, exist_ok=True)
	for res, normalization in resonances.items():
		if res in grids:
			output_path = output_directory / f"{res}.pineappl.lz4"
			print(f"Converting {res} and writing it to {output_path}")
			convert_grid_to_pineappl(grids.pop(res), normalization).write_lz4(output_path)

def dump_grid_directory_to_pineappl(grid_directory:Path, output_directory:Path, resonances:dict[str, float]={"D0": 1.0, "Dp": 1.0, "Ds": 1.0, "LambdaC": 1.14}, pattern:str = "*.csv.gz") -> None:
	"""
	Converts all grid files in the given directory to PineAPPL grids, one grid file at a time

	Parameters
	----------
	grid_directory:Path
		Path pointing to the directory containing the grid files
	output_directory:Path
		Path pointing to the directory in which the PineAPPL grids will be written
	resonances:dict[str, float] [default: {"D0": 1.0, "Dp": 1.0, "Ds": 1.0, "LambdaC": 1.14}]
		Dictionary defining the resonance particle names (which must match the ones in the input grids) and their corresponding normalization constants,
		which will multiply all appropriate values in the output grids
	pattern:str	[default: "*.csv.gz"]
		Pattern string used for globbing the files in `grid_directory`

	Notes
	-----
	The PineAPPL grids of each grid file are written by `convert_gridfile_to_pineappl` to the subdirectory '<output_directory>/<grid file>', where <grid file> is the name of the grid file without the '.csv.gz' suffix
	"""

	print(f"Searching directory \n  {grid_directory}\nwith pattern '{pattern}'... ", end="")
	input_paths = sorted(grid_directory.glob(pattern))
	print(f"found {len(input_paths)} grid files")

	if not input_paths:
		raise FileNotFoundError(f"No grid files matching '{pattern}' in {grid_directory}")

	for input_path in input_paths:
		convert_gridfile_to_pineappl(input_path, output_directory / input_path.name.removesuffix(".csv.gz"), resonances)
	print("All done")