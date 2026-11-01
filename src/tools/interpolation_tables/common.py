from __future__ import annotations

import itertools
from dataclasses import dataclass
import pandas as pd
import numpy as np

### Utils

type float_array = np.ndarray[tuple[int], np.dtype[np.floating]]

@dataclass
class interpolation_grid:
	df:pd.DataFrame
	acc:accelerator

@dataclass
class accelerator:
	unique_xs:float_array
	unique_ys:float_array
	unique_Es:float_array
	xs:float_array
	ys:float_array
	Es:float_array
	Q2s:float_array
	renormalization_scales:dict
	factorization_scales:dict
	fragmentation_scales:dict
	initial_flavors:float_array
	fragmentation_flavors:float_array
	x_knot_indices:float_array
	z_knot_indices:float_array
	x_knots:float_array
	z_knots:float_array
	unique_resonances:list[str]

def make_accelerator(df:pd.DataFrame) -> accelerator:
	representative_slice = get_representative_slice(df)

	grouping = representative_slice.groupby(["y", "E"])

	def make_dictionary(col:str) -> dict:
		keys = grouping.groups.keys()
		values = grouping.apply(lambda g: g[col].values, include_groups=False)
		return dict(zip(keys, values))

	return accelerator(
		get_unique_values(representative_slice, "x"),
		get_unique_values(representative_slice, "y"),
		get_unique_values(representative_slice, "E"),
		slice_to_array(representative_slice, "x"),
		slice_to_array(representative_slice, "y"),
		slice_to_array(representative_slice, "E"),
		slice_to_array(representative_slice, "Q2"),
		make_dictionary("ren"),
		make_dictionary("fact"),
		make_dictionary("frag"),
		get_unique_values(df, "if"),
		get_unique_values(df, "ff"),
		get_unique_values(df, "ik"),
		get_unique_values(df, "jk"),
		get_unique_values(df, "xk"),
		get_unique_values(df, "zk"),
		get_unique_values(df, "res")
	)

def slice_to_array(slice:pd.DataFrame, column:str) -> float_array:
	return slice[column].to_numpy()

def get_unique_values(df:pd.DataFrame, col:str) -> float_array:
	return df[col].unique()

def get_representative_slice(df:pd.DataFrame) -> pd.DataFrame:
	i_f = get_unique_values(df, "if")[0]
	f_f = get_unique_values(df, "ff")[0]
	ik = get_unique_values(df, "ik")[0]
	jk = get_unique_values(df, "jk")[0]

	representative_slice = df[np.logical_and.reduce([df["if"] == i_f, df["ff"] == f_f, df["ik"] == ik, df["jk"] == jk])]
	# A grid that is not split by resonance has each kinematical point once for every resonance, so only the first resonance is kept
	return representative_slice[representative_slice["res"] == representative_slice["res"].iloc[0]]

### Dataframe concatination

def import_gridfile(filename:str, compression={"method": "gzip", "compresslevel": 1}, fill_missing_channels:bool=True, split_by_resonance:bool=True) -> list[interpolation_grid] | interpolation_grid:
	"""
	With `split_by_resonance`, the grid file is returned as a list of grids, one for each resonance. The convolutions in `convolution.py`
	and the PineAPPL conversion need this. Otherwise the whole grid file is returned as a single grid.
	"""
	df = pd.read_csv(filename, compression=compression, comment="#").sort_values(["y", "E", "if", "ff", "x", "ik", "jk"])
	if not split_by_resonance:
		return import_grid(df, fill_missing_channels)
	resonances = get_unique_values(df, "res")
	dfs_by_resonance = [df[df["res"] == res] for res in resonances]
	# The slices are copies, so the full dataframe is no longer needed
	del df
	# Each slice is released as soon as its grid is built
	grids_by_resonance = []
	while dfs_by_resonance:
		grids_by_resonance.append(import_grid(dfs_by_resonance.pop(0), fill_missing_channels))
	return grids_by_resonance

def import_grid(df:pd.DataFrame, fill_missing_channels:bool=True) -> interpolation_grid:
	"""
	With `fill_missing_channels`, every combination of initial and fragmentation flavors missing from the grid is added with zero values.
	The convolutions in `convolution.py` need this, the PineAPPL conversion does not.
	"""
	acc = make_accelerator(df)

	if not fill_missing_channels:
		return interpolation_grid(df, acc)

	all_flavor_combinations = itertools.product(acc.initial_flavors, acc.fragmentation_flavors)
	representative = df[np.logical_and(df["if"] == acc.initial_flavors[0], df["ff"] == acc.fragmentation_flavors[0])]

	def build_new_df(df):
		def build_dummy_df(i_f, f_f):
			dummy_df = representative.copy()
			dummy_df["LO"] = 0.0
			dummy_df["NLO_correction"] = 0.0
			dummy_df["NNLO_correction"] = 0.0
			dummy_df["if"] = i_f
			dummy_df["ff"] = f_f
			return dummy_df

		dummies = [build_dummy_df(i_f, f_f) for i_f, f_f in all_flavor_combinations if df[np.logical_and(df["if"] == i_f, df["ff"] == f_f)].empty]
		new_df = pd.concat([df, *dummies]).sort_values(["y", "E", "if", "ff", "x", "ik", "jk"])
		return new_df

	new_df = build_new_df(df)
	new_acc = make_accelerator(new_df)

	return interpolation_grid(new_df, new_acc)