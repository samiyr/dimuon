import subprocess
import itertools
import numpy as np
import pandas as pd
from io import StringIO
from .arg_modes import Observable
from .arg_modes import RunMode
from .arg_modes import ResonanceParticle
from . import arg_parser

from typing import Optional

def dimuon_cross_section(arg_list):
	"""
	Parses the list of arguments, calls the external dimuon program in direct mode to perform the actual cross-section computation,
	and parses the result into a list of output values. By default, the dimuon program is the one installed alongside
	the pool controller (see --cmd).
	"""
	args 											= arg_list[0]
	x:float 										= arg_list[1]
	y:float 										= arg_list[2]
	E:float 										= arg_list[3]
	pdf_member:int 									= arg_list[4]
	neutrino:bool 									= arg_list[5]
	kinematic_folder:str 							= arg_list[6]
	charm_mass:float 								= arg_list[7]
	(renormalization, factorization, fragmentation) = arg_list[8]
	parametrization:Optional[tuple[float]] 			= arg_list[9]

	command_arguments = [
		args.cmd, 
		"--observable", str(args.observable), 
		"--pdf", args.pdf,
		"--pdf-member", str(pdf_member),
		"--x", str(x),
		"--y", str(y),
		"--E", str(E),
		"--renormalization", str(renormalization),
		"--factorization", str(factorization),
		"--fragmentation", str(fragmentation),
		"--charm-flavors", args.charm_flavors,
		"--integration-method", str(args.integration_method),
		"--cuba-min-evals", str(args.cuba_min_evals),
		"--cuba-max-rel-error-lo", str(args.cuba_max_rel_error_lo),
		"--cuba-max-rel-error-nlo", str(args.cuba_max_rel_error_nlo),
		"--cuba-max-rel-error-nnlo", str(args.cuba_max_rel_error_nnlo),
		"--cuba-suave-subd-evals", str(args.cuba_suave_subd_evals),
		"--cuba-suave-min-samples", str(args.cuba_suave_min_samples),
		"--cuba-suave-flatness", str(args.cuba_suave_flatness),
		"--cuba-vegas-start-evals", str(args.cuba_vegas_start_evals),
		"--cuba-vegas-increase-evals", str(args.cuba_vegas_increase_evals)
	] + arg_parser.shared_command_arguments(args, RunMode.direct, neutrino)

	if charm_mass != None:
		# If an explicit charm mass is requested, pass it on. Otherwise, don't pass a charm mass argument at all,
		# which signals to the CC program that it should get the charm mass from the PDF.
		command_arguments += ["--charm-mass", np.format_float_positional(charm_mass)]

	if args.isospin:
		command_arguments += ["--isospin"]
		command_arguments += ["--Z", str(args.Z)]
		command_arguments += ["--A", str(args.A)]

	if parametrization != None:
		(N, alpha, beta, gamma) = parametrization
		command_arguments += ["--decay-N", str(N)]
		command_arguments += ["--decay-alpha", str(alpha)]
		command_arguments += ["--decay-beta", str(beta)]
		command_arguments += ["--decay-gamma", str(gamma)]

	# Launch the CC program
	result = subprocess.run(
		command_arguments, 
		stdout=subprocess.PIPE, # for obtaining the path of the output file
		check=True
	)

	# Split the output into an array of column values.
	# They are on the line after the first "===", followed by the channel decomposition (for both SIDIS and DIS).
	output_list = arg_parser.read_output_file(result.stdout)
	result_line = output_list[output_list.index("===") + 1]
	values = result_line.split(", ")

	# Include some internal data not returned by the CC program, but that is needed for bookkeeping later
	values += [kinematic_folder, charm_mass]
	values += [(renormalization, factorization, fragmentation)]
	values += [parametrization]
	return values

def list_to_string(l:list) -> str:
	return ','.join(map(str, l))

def interpolation_table(arg_list):
	"""
	Parses the list of arguments, calls the external dimuon program in kernel mode to compute the rows of the interpolation table
	at one kinematic point for one resonance particle, and parses the result into a dataframe with one row per pair of kernels and physical channel.
	"""
	args						= arg_list[0]
	x:float						= arg_list[1]
	y:float						= arg_list[2]
	E:float						= arg_list[3]
	resonance:ResonanceParticle	= arg_list[4]
	pdf_knots:np.ndarray		= arg_list[5]
	ff_knots:np.ndarray			= arg_list[6]

	command_arguments = [
		args.cmd,
		"--observable", str(Observable.sidis),
		"--resonance-particle", str(resonance),
		"--pdf-knots", list_to_string(pdf_knots),
		"--ff-knots", list_to_string(ff_knots),
		"--x", str(x),
		"--y", str(y),
		"--E", str(E),
		# The CC program has no PDF in the kernel mode, so the charm mass has to be given
		"--charm-mass", np.format_float_positional(args.charm_mass)
	] + arg_parser.shared_command_arguments(args, RunMode.kernel, args.neutrino)

	# Launch the CC program
	result = subprocess.run(
		command_arguments, 
		stdout=subprocess.PIPE, # for obtaining the path of the output file
		check=True
	)

	# The output file has the kinematics and the scales, and then one row per pair of kernels and physical channel, separated by "===" lines
	output_list = arg_parser.read_output_file(result.stdout)

	result_line, channel_lines = [list(y) for x, y in itertools.groupby(output_list, lambda z: z == "===") if not x]

	results = result_line[0].split(",")

	df = pd.read_csv(StringIO("\n".join(channel_lines)), delimiter=",", names=["ik", "jk", "if", "ff", "LO", "NLO_correction", "NNLO_correction"])

	df["res"] = str(resonance)
	df["xk"] = pdf_knots[df["ik"].to_numpy()]
	df["zk"] = ff_knots[df["jk"].to_numpy()]
	df["x"] = x
	df["y"] = y
	df["E"] = E
	df["Q2"] = float(results[3])
	df["ren"] = float(results[4])
	df["fact"] = float(results[5])
	df["frag"] = float(results[6])

	return df[["if", "ff", "LO", "NLO_correction", "NNLO_correction", "res", "ik", "xk", "jk", "zk", "x", "y", "E", "Q2", "ren", "fact", "frag"]]

def pion_production(arg_list):
	"""
	Parses the list of arguments, calls the external pionprod-cp program for one bin in z and Q2, and parses the
	output into a list of output values. The pionprod-cp program must be located in the directory from where
	the pion_cross_sections.py program is launched.
	"""
	z_min, z_max, Q2_min, Q2_max, Wsign = arg_list

	command_arguments = [
		"./pionprod-cp",
		str(z_min), str(z_max),
		str(Q2_min), str(Q2_max),
		str(Wsign)
	]

	result = subprocess.run(
		command_arguments, 
		stdout=subprocess.PIPE, # for obtaining the output
		check=True
	)

	# Split the output string into an array of column values
	output_list = result.stdout.decode("utf-8").splitlines()
	values = output_list[-1].split(", ")

	return values