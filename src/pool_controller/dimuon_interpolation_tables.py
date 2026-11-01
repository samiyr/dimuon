import numpy as np
from scipy.special import lambertw
import itertools
import functools
import operator
from utils import constants
import time
import datetime
from utils import launch
from utils import arg_parser
from utils import metadata
from utils import pool
import sys
import os
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import gzip
from pathlib import Path
from utils.arg_modes import RunMode

def lambertgrid(x_min, x_max, n_points, reflected=False):
	def direct_relation(x):
		if reflected:
			return 5 * x - np.log(1 - x)
		else:
			return 5 * (1 - x) - np.log(x)

	def inverse_relation(y):
		if reflected:
			return np.real((1 / 5) * (5 - lambertw(5 * np.exp(5 - y))))
		else:
			return np.real(1 / 5 * lambertw(5 * np.exp(5 - y)))
			
	y_min = direct_relation(x_min)
	y_max = direct_relation(x_max)

	return np.array([inverse_relation(y) for y in np.linspace(y_min, y_max, n_points)])

if __name__ == '__main__':
	parser = arg_parser.create_arg_parser(RunMode.kernel)

	args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
	
	if args.neutrino and args.antineutrino:
		raise ValueError("The flags --neutrino and --antineutrino cannot be used simultaneously with dimuon_interpolation_tables.py")
	elif not args.neutrino and not args.antineutrino:
		raise ValueError("One of the flags --neutrino or --antineutrino has to be included with dimuon_interpolation_tables.py")

	if args.nutev and args.ccfr:
		raise ValueError("The flags --nutev and --ccfr cannot be used simultaneously with dimuon_interpolation_tables.py")
	elif not args.nutev and not args.ccfr:
		raise ValueError("One of the flags --nutev or --ccfr has to be included with dimuon_interpolation_tables.py")

	# The CC program computes NNLO if there are NNLO channels, and then it needs the thresholds, since it has no PDF to take them from
	nnlo_channels = any(argument.startswith("--nnlo-") for argument in constants.get_channels(args.neutrino, args))
	if nnlo_channels and None in [args.charm_threshold, args.bottom_threshold, args.top_threshold]:
		raise ValueError(
			"The flags --charm-threshold, --bottom-threshold and --top-threshold all have to be given with dimuon_interpolation_tables.py if there are NNLO channels, "
			"since they determine the number of active flavors in the NNLO coefficients. "
			"They should match the quark masses of the PDF set that the grids are convolved with."
		)

	pdf_knots = lambertgrid(1e-2, 9.999e-1, args.pdf_knot_count, reflected=False)
	ff_knots = lambertgrid(1e-2, 9.9989e-1, args.ff_knot_count, reflected=True)

	# Kinematical arguments based on passed args
	if args.nutev and args.neutrino:
		xs = constants.nutev_neutrino_x
		ys = constants.nutev_neutrino_y
		Es = constants.nutev_neutrino_E
	elif args.nutev and args.antineutrino:
		xs = constants.nutev_antineutrino_x
		ys = constants.nutev_antineutrino_y
		Es = constants.nutev_antineutrino_E
	elif args.ccfr and args.neutrino:
		xs = constants.ccfr_neutrino_x
		ys = constants.ccfr_neutrino_y
		Es = constants.ccfr_neutrino_E
	elif args.ccfr and args.antineutrino:
		xs = constants.ccfr_antineutrino_x
		ys = constants.ccfr_antineutrino_y
		Es = constants.ccfr_antineutrino_E

	if not args.experimental_x:
		xs = constants.x_bins

	if args.x:
		xs = args.x

	if args.y:
		ys = args.y
	
	if args.E:
		Es = args.E
	
	# Generate run args
	iters = [
		[args],
		xs, ys, Es,
		args.resonance_particles,
		[pdf_knots], [ff_knots]
	]
	run_args = itertools.product(*iters)

	print(f"Constructed {functools.reduce(operator.mul, map(len, iters), 1)} computations")

	start = time.time()
	dfs = pool.run(launch.interpolation_table, run_args, mpi=args.mpi)

	print("Creating dataframe")
	df = pd.concat(dfs, ignore_index=True).sort_values(["y", "E", "if", "ff", "x", "ik", "jk"])

	print("Dataframe shape:", df.shape)

	def export(df:pd.DataFrame, 
		filename:str, dir:Path
	):
		"""
		Export the supplied dataframe as a gzip-compressed CSV file, named filename with ".gz" appended, to a given directory. 
		Generates the directory structure as needed.
		"""

		dir.mkdir(parents=True, exist_ok=True)
		file_path = dir / f"{filename}.gz"

		# The metadata comments and the dataframe are written to the same gzip stream, which compresses them on the fly

		def write_df(file):
			print(f"Writing dataframe to file {file_path}")

			# The index only has the row numbers from before the sort, so it is not written
			pa_df = pa.Table.from_pandas(df, preserve_index=False)
			csv.write_csv(pa_df, file)

		with gzip.open(file_path, 'wt', compresslevel=1) as file:
			metadata.write(file, args)

		with gzip.open(file_path, 'ab', compresslevel=1) as file:
			write_df(file)

		print("Dataframe written to file")

	# Read the environment variable where the output directory structure should be generated in
	global_dir = Path(os.environ.get("DIMUON_DATA_DIR", ""))

	data_folder = "kernel-mode"
	scheme_folder = args.dir if args.dir != None else ""

	export(df, args.grid_output, global_dir / data_folder / scheme_folder)

	delta_seconds = time.time() - start
	delta = datetime.timedelta(seconds=int(delta_seconds))
	print("Done, took", str(delta))