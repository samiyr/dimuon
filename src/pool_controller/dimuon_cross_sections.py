from utils import launch
from utils import scale_variations
from utils import constants
from utils import arg_generator
from utils import metadata
from utils import pool

from utils import arg_parser

from utils.arg_modes import Mode, RunMode

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import time
import datetime
from typing import Union, Optional

def generate_run_args(args):
	# Generate run arguments
	run_args = []

	if args.y and args.E:
		def custom_kinematics(neutrino:bool):
			return arg_generator.differential(
				args,
				args.x if args.x else constants.x_bins, args.y, args.E, neutrino, args.kinematic_folder
			)
			
		if args.neutrino:
			run_args += custom_kinematics(True)
		if args.antineutrino:
			run_args += custom_kinematics(False)

	if (args.nutev and args.neutrino):
		run_args += arg_generator.differential(
			args,
			constants.x_bins if not args.experimental_x else constants.nutev_neutrino_x, 
			constants.nutev_neutrino_y, constants.nutev_neutrino_E, True, "nutev_neutrino"
		)
	if (args.nutev and args.antineutrino):
		run_args += arg_generator.differential(
			args,
			constants.x_bins if not args.experimental_x else constants.nutev_antineutrino_x, 
			constants.nutev_antineutrino_y, constants.nutev_antineutrino_E, False, "nutev_antineutrino"
		)
	if (args.ccfr and args.neutrino):
		run_args += arg_generator.differential(
			args,
			constants.x_bins if not args.experimental_x else constants.ccfr_neutrino_x, 
			constants.ccfr_neutrino_y, constants.ccfr_neutrino_E, True, "ccfr_neutrino"
		)
	if (args.ccfr and args.antineutrino):
		run_args += arg_generator.differential(
			args,
			constants.x_bins if not args.experimental_x else constants.ccfr_antineutrino_x, 
			constants.ccfr_antineutrino_y, constants.ccfr_antineutrino_E, False, "ccfr_antineutrino"
		)

	return run_args

def handle_output(outputs, args):
	# By default, no scale variations
	scales = [(1.0, 1.0, 1.0)]
	# If the run mode is scales, then include all scale variations
	if args.mode == Mode.scales7p or args.mode == Mode.scales17p:
		scales = list(scale_variations.generate_standard(args))	
		
	parametrizations = [None]
	if args.mode == Mode.decay_function_errors and args.decay_parametrizations:
		parametrization_list = np.loadtxt(args.decay_parametrizations, delimiter=",")
		parametrizations = [tuple(p) for p in parametrization_list]

	# Create a dataframe from output list
	print("Creating dataframe")
	df = pd.DataFrame(outputs, columns=constants.get_df_columns(args.observable))

	print("Dataframe shape:", df.shape)

	# Read the environment variable where the output directory structure should be generated in
	global_dir = Path(os.environ.get("DIMUON_DATA_DIR", ""))

	data_folder = constants.get_data_folder(args)
	scheme_folder = args.dir if args.dir != None else ""

	output_columns = constants.get_output_columns(args)

	def export(df:pd.DataFrame, 
			filename:str, dir:Path, 
			pdf:str, member:int, 
			charm_mass:Union[float, None] = None, 
			renormalization:float = 1.0, factorization:float = 1.0, fragmentation:float = 1.0,
			N:Optional[float] = None, alpha:Optional[float] = None, beta:Optional[float] = None, gamma:Optional[float] = None
		):
		"""
		Export the supplied dataframe as a CSV file to a given directory. 
		Generates the directory structure as needed.
		"""

		dir.mkdir(parents=True, exist_ok=True)
		file_path = dir / filename

		with open(file_path, 'w') as file:
			metadata.write(file, args)

			file.write(f"# charm mass: {charm_mass}\n")
			file.write(f"# pdf set: {pdf}\n")
			file.write(f"# pdf member: {member}\n")
			file.write(f"# renormalization: {renormalization}\n")
			file.write(f"# factorization: {factorization}\n")
			file.write(f"# fragmentation: {fragmentation}\n")

			if N != None and alpha != None and beta != None and gamma != None:
				file.write(f"# N: {N}\n")
				file.write(f"# alpha: {alpha}\n")
				file.write(f"# beta: {beta}\n")
				file.write(f"# gamma: {gamma}\n")

		df.to_csv(dir / filename, index=False, header=True, columns=output_columns, mode="a", float_format="%.12f")

	# Split the single output dataframe into multiple ones based on PDF set, category (i.e. kinematic folder), and PDF member
	print("Splicing output data")
	for pdf, df_pdf in df.groupby("pdf_set"):
		for category, df_category in df_pdf.groupby("category"):
			for member, df_member in df_category.groupby("pdf_member"):
				member_filename = f"{member.zfill(4)}.csv"

				if args.mode == Mode.normal:
					export(
						df_member, member_filename, 
						global_dir / data_folder / pdf / scheme_folder / category, 
						pdf=pdf, member=member
					)
				elif args.mode == Mode.masses:
					for mass, df_mass in df_member.groupby("charm_mass"):
						mass_folder = np.format_float_positional(mass).replace(".", "_")
						export(
							df_mass, member_filename, 
							global_dir / data_folder / pdf / scheme_folder / mass_folder / category,
							pdf=pdf, member=member, charm_mass=mass
						)
				elif args.mode == Mode.scales7p or args.mode == Mode.scales17p:
					for (renormalization, factorization, fragmentation), df_scale in df_member.groupby("scales"):
						scale_index = scales.index((renormalization, factorization, fragmentation))
						scale_filename = "base_scale.csv" if renormalization == 1.0 and factorization == 1.0 and fragmentation == 1.0 else f"scale_{scale_index}.csv"
						export(
							df_scale, scale_filename, global_dir / data_folder / pdf / scheme_folder / category,
							pdf=pdf, member=member, renormalization=renormalization, factorization=factorization, fragmentation=fragmentation
						)
				elif args.mode == Mode.decay_function_errors:
					for (N, alpha, beta, gamma), df_parametrization in df_member.groupby("parametrization"):
						parametrization_index = parametrizations.index((N, alpha, beta, gamma))
						parametrization_filename = f"{str(parametrization_index).zfill(4)}.csv"
						export(
							df_parametrization, parametrization_filename, global_dir / data_folder / pdf / scheme_folder / category,
							pdf=pdf, member=member, N=N, alpha=alpha, beta=beta, gamma=gamma
						)

if __name__ == '__main__':
	args = arg_parser.create_arg_parser(RunMode.direct).parse_args(args=None if sys.argv[1:] else ['--help'])
	print("Parsed arguments:", args)

	if not args.neutrino and not args.antineutrino:
		raise ValueError("One of the flags --neutrino or --antineutrino has to be included with dimuon_cross_sections.py")

	if bool(args.y) != bool(args.E):
		raise ValueError("The flags --y and --E have to be included together with dimuon_cross_sections.py")
	elif args.y and not args.kinematic_folder:
		raise ValueError("The flag --kinematic-folder has to be included with --y and --E")

	if not args.nutev and not args.ccfr and not args.y:
		raise ValueError("One of the flags --nutev, --ccfr, or --y and --E has to be included with dimuon_cross_sections.py")

	run_args = generate_run_args(args)
	print("Generated run arguments, count", len(run_args))

	start = time.time()
	outputs = pool.run(launch.dimuon_cross_section, run_args, mpi=args.mpi)

	handle_output(outputs, args)

	delta_seconds = time.time() - start
	delta = datetime.timedelta(seconds=int(delta_seconds))
	print("Done, took", str(delta))