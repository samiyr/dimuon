from .arg_modes import Mode
from .arg_modes import Observable
from .arg_modes import RunMode
from .arg_modes import ScaleDependence
from .arg_modes import Channels
from .arg_modes import CoefficientImplementation
from .arg_modes import IntegrationMethod
from .arg_modes import ResonanceParticle

from . import constants

import argparse
import configargparse
import numpy as np
import os
from pathlib import Path

def create_arg_parser(run_mode:RunMode):
	"""
	Creates the argument parser of dimuon_cross_sections.py (the direct mode) or of dimuon_interpolation_tables.py (the kernel mode), with the arguments of that run mode.
	The arguments that are passed on to the CC program have the same names there.

	Args:
		run_mode: Run mode of the CC program
	"""
	# Abbreviated arguments are not allowed, since they would take e.g. --pdf as --pdf-knot-count in dimuon_interpolation_tables.py
	parser = configargparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		config_file_parser_class=configargparse.YAMLConfigFileParser,
		allow_abbrev=False
	)

	float_list_parser = lambda t: [float(s.strip()) for s in t.split(",")] if t else []

	# The default command is the installed dimuon executable.
	# This file is installed to <installdir>/opt/dimuon/pool_controller/utils, and the executable to <installdir>/bin.
	installed_cmd = Path(__file__).resolve().parents[4] / "bin" / "dimuon"
	parser.add_argument("--cmd", type=str, default=str(installed_cmd), help="command to run")

	# Techical options
	parser.add_argument("--config", is_config_file=True, help="run card file")
	parser.add_argument("--mpi", action="store_true", default=False, help="use MPI (requires mpi4py package)")
	parser.add_argument("--version", type=str, default="0", help="specify version of the run")

	# Kinematics
	parser.add_argument("--nutev", action="store_true", default=False, help="include NuTeV kinematics")
	parser.add_argument("--ccfr", action="store_true", default=False, help="include CCFR kinematics")
	parser.add_argument("--experimental-x", action="store_true", default=False, help="use only experimental Bjorken-x values")
	parser.add_argument("--x", type=float_list_parser, default="", help="optional list of comma-separated values to use as Bjorken-x values instead of defaults")
	parser.add_argument("--y", type=float_list_parser, default="", help="optional list of comma-separated values to use as inelasticities instead of defaults")
	parser.add_argument("--E", type=float_list_parser, default="", help="optional list of comma-separated values to use as beam energies instead of defaults")
	if run_mode == RunMode.direct:
		parser.add_argument("--kinematic-folder", type=str, help="optional name of kinematic folder, used in conjunction with --y and --E")

	# Process selection
	parser.add_argument("--neutrino", action="store_true", default=False, help="include neutrino scattering")
	parser.add_argument("--antineutrino", action="store_true", default=False, help="include antineutrino scattering")

	# Channel options. The CC program computes up to the highest perturbative order that has channels.
	channel_group = parser.add_mutually_exclusive_group(required=True)
	channel_group.add_argument("--channels", type=Channels, choices=list(Channels), help="channel preset, computed up to the highest perturbative order that has channels")
	channel_group.add_argument("--channel-list", action="append", help="explicit list of channel configurations, computed up to the highest perturbative order that has channels")
	if run_mode == RunMode.direct:
		parser.add_argument("--charm-flavors", type=str, default="1,3,5,4", help="list of comma-separated flavors for inclusive charm production")
	parser.add_argument("--flavors", type=str, default="1,3,5,2,4", help="list of comma-separated flavors for semi-inclusive dimuon production")
	parser.add_argument("--initial-flavors", type=str, required=False, help="list of comma-separated flavors used for filtering SIDIS channels based on the initial-state flavors")
	parser.add_argument("--fragmentation-flavors", type=str, required=False, help="list of comma-separated flavors used for filtering SIDIS channels based on the final-state fragmenting flavors")

	# Observable options. Gridding is for the SIDIS observable only.
	if run_mode == RunMode.direct:
		parser.add_argument("--observable", type=Observable, choices=list(Observable), required=True, help="observable to compute")
	parser.add_argument("--coefficient-implementation-lo", type=CoefficientImplementation, choices=list(CoefficientImplementation), default=CoefficientImplementation.hw, help="coefficient implementation at LO")
	parser.add_argument("--coefficient-implementation-nlo", type=CoefficientImplementation, choices=list(CoefficientImplementation), default=CoefficientImplementation.hw, help="coefficient implementation at NLO")
	parser.add_argument("--coefficient-implementation-nnlo", type=CoefficientImplementation, choices=list(CoefficientImplementation), default=CoefficientImplementation.gc, help="coefficient implementation at NNLO (only gc)")

	# Scale dependence options
	parser.add_argument("--scale", type=ScaleDependence, choices=list(ScaleDependence), default=ScaleDependence.Q2mc2, help="scale dependence")
	parser.add_argument("--renormalization-unbounded", action="store_true", default=False, help="if enabled, renormalization scale will not be bounded from below by the charm mass")
	parser.add_argument("--factorization-unbounded", action="store_true", default=False, help="if enabled, factorization scale will not be bounded from below by the charm mass")
	parser.add_argument("--fragmentation-unbounded", action="store_true", default=False, help="if enabled, fragmentation scale will not be bounded from below by the charm mass")

	# Mass options. Gridding has no PDF to take the charm mass from, so there it is required.
	if run_mode == RunMode.direct:
		parser.add_argument("--masses", type=float_list_parser, default="0.01,1.3,1.5", help="list of comma-separated values used as charm masses in --mode=masses")
		parser.add_argument("--charm-mass", type=float, default=None, help="override the charm mass of the PDF in modes other than --mode=masses")
	else:
		parser.add_argument("--charm-mass", type=float, required=True, help="charm mass")

	# Flavor thresholds, used for the number of active flavors in the NNLO coefficients.
	# These are separate from the quark masses above, which set the final-state masses of the channels.
	# By default they are the quark masses of the PDF. Gridding has no PDF, so there they are required if there are NNLO channels.
	threshold_help = "override the {0} flavor threshold (default: the {0} mass of the PDF)" if run_mode == RunMode.direct \
		else "{0} flavor threshold, required if there are NNLO channels: the {0} mass of the PDF sets that the grids are convolved with"
	parser.add_argument("--charm-threshold", type=float, default=None, help=threshold_help.format("charm"))
	parser.add_argument("--bottom-threshold", type=float, default=None, help=threshold_help.format("bottom"))
	parser.add_argument("--top-threshold", type=float, default=None, help=threshold_help.format("top"))

	if run_mode == RunMode.direct:
		# PDF
		parser.add_argument("--pdf", type=str, required=True, help="name of PDF set")
		parser.add_argument("--pdf-start", type=int, default=0, help="first PDF member to include")
		parser.add_argument("--pdf-end", type=int, default=0, help="last PDF member to include")
		parser.add_argument("--isospin", action="store_true", default=False, help="enable isospin corrections")
		parser.add_argument("--Z", type=float, default=1.0, help="atomic number used when isospin corrections are enabled")
		parser.add_argument("--A", type=float, default=1.0, help="mass number used when isospin corrections are enabled")

		# Variation options
		parser.add_argument("--mode", type=Mode, choices=list(Mode), default=Mode.normal, help="variations to compute")
		parser.add_argument("--dir", type=str, default=None, help="if given, this additional directory name will be introduced to the directory structure when saving output")

		# Decay function parametrization
		parser.add_argument("--decay-parametrizations", type=str, default=None, help="path to decay function parametrization file")
		parser.add_argument("--decay-start", type=int, default=0, help="first decay parametrization to include")
		parser.add_argument("--decay-end", type=int, default=0, help="last decay parametrization to include")

	# Integration options. The CC program integrates the direct mode with --integration-method, and the kernel mode (gridding) with plain quasi-Monte Carlo.
	# The default evaluations are a power of two, which suits the Sobol points of plain quasi-Monte Carlo.
	parser.add_argument("--cuba-max-evals", type=int, default=2**20, help="maximum number of evaluations, or the number of quasi-Monte Carlo points")
	if run_mode == RunMode.direct:
		parser.add_argument("--integration-method", type=IntegrationMethod, choices=list(IntegrationMethod), default=IntegrationMethod.suave, help="integration method")
		# The Suave and Vegas iteration sizes are 100k, so with the default evaluations they run 11 iterations
		parser.add_argument("--cuba-min-evals", type=int, default=0)
		parser.add_argument("--cuba-max-rel-error-lo", type=float, default=1e-4)
		parser.add_argument("--cuba-max-rel-error-nlo", type=float, default=1e-3)
		parser.add_argument("--cuba-max-rel-error-nnlo", type=float, default=1e-3)
		parser.add_argument("--cuba-suave-subd-evals", type=int, default=100_000)
		parser.add_argument("--cuba-suave-min-samples", type=int, default=25)
		parser.add_argument("--cuba-suave-flatness", type=float, default=50.0)
		parser.add_argument("--cuba-vegas-start-evals", type=int, default=100_000, help="evaluations per Vegas iteration")
		parser.add_argument("--cuba-vegas-increase-evals", type=int, default=0, help="increase in evaluations per Vegas iteration")

	if run_mode == RunMode.kernel:
		# Gridding options
		parser.add_argument("--grid-output", type=str, required=True, help="name of the grid output file")
		parser.add_argument("--dir", type=str, default=None, help="if given, this additional directory name will be introduced to the directory structure when saving the grid output file")
		parser.add_argument("--resonance-particles", type=lambda t: [ResonanceParticle(s.strip()) for s in t.split(",")] if t else [], default="D0,Dp,Ds,LambdaC", help="resonance particles for gridding")
		parser.add_argument("--pdf-knot-count", type=int, required=True, help="number of PDF interpolation knots")
		parser.add_argument("--ff-knot-count", type=int, required=True, help="number of FF interpolation knots")

	return parser

def shared_command_arguments(args, run_mode:RunMode, neutrino:bool) -> list[str]:
	"""
	Returns the arguments of the CC program that dimuon_cross_sections.py and dimuon_interpolation_tables.py pass in the same way, from the arguments of both run modes.

	Args:
		args: Arguments parsed by the parser of create_arg_parser
		run_mode: Run mode of the CC program
		neutrino: If true, do neutrino scattering. Otherwise, do antineutrino scattering
	"""
	command_arguments = [
		"--run-mode", str(run_mode),
		"--output-mode", "file",
		"--neutrino" if neutrino else "--antineutrino",
		"--scale", str(args.scale),
		"--coefficient-implementation-lo", str(args.coefficient_implementation_lo),
		"--coefficient-implementation-nlo", str(args.coefficient_implementation_nlo),
		"--coefficient-implementation-nnlo", str(args.coefficient_implementation_nnlo),
		"--flavors", args.flavors,
		"--cuba-max-evals", str(args.cuba_max_evals)
	]

	# Thresholds of the heavy quarks, which determine the number of active flavors in the NNLO coefficients.
	# If they are not given, the CC program takes them from the PDF.
	for threshold, flag in [
		(args.charm_threshold, "--charm-threshold"),
		(args.bottom_threshold, "--bottom-threshold"),
		(args.top_threshold, "--top-threshold")
	]:
		if threshold is not None:
			command_arguments += [flag, np.format_float_positional(threshold)]

	if args.renormalization_unbounded:
		command_arguments += ["--renormalization-unbounded"]
	if args.factorization_unbounded:
		command_arguments += ["--factorization-unbounded"]
	if args.fragmentation_unbounded:
		command_arguments += ["--fragmentation-unbounded"]

	if args.initial_flavors:
		command_arguments += ["--initial-flavors", args.initial_flavors]
	if args.fragmentation_flavors:
		command_arguments += ["--fragmentation-flavors", args.fragmentation_flavors]

	command_arguments += constants.get_channels(neutrino, args)

	return command_arguments

def read_output_file(stdout:bytes) -> list[str]:
	"""
	Returns the lines of the output file of the CC program, and deletes the file.
	With --output-mode file (see shared_command_arguments), the CC program writes its output to a file
	and prints the path of the file on the last line of stdout.

	Args:
		stdout: Standard output of the CC program
	"""
	output_file = stdout.decode("utf-8").splitlines()[-1]

	with open(output_file) as file:
		output_lines = file.read().splitlines()
	os.remove(output_file)

	return output_lines