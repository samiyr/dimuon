from .arg_modes import Mode

from typing import Optional
import itertools
from . import scale_variations

import numpy as np

def full_differential(
		args,
		xs: list[float], ys: list[float], Es: list[float], 
		start: int, end: int,
		neutrino: bool,
		kinematic_folder: str,
		charm_masses: list[Optional[float]], scales: list[tuple[float]],
		parametrizations: list[Optional[tuple[float]]]
	):
	"""
	Generate an iterable of all argument combinations.

	Args:
		args: Instance of argparse
		xs: List of Bjorken-x values
		ys: List of inelasticity (y) values
		Es: List of beam energy (E) values
		start: First PDF member to include
		end: last PDF member to include
		neutrino: If true, do neutrino scattering. Otherwise, do antineutrino scattering
		kinematic_folder: Name of the kinematic folder, used for generating the output directory structure
		charm_masses: List of charm mass values. A value of None indicates that the mass should be read from the PDF
		scales: List of scale multipliers, organized into a tuple (renormalization, factorization, fragmentation)
	"""
	return itertools.product(
		[args],
		xs, ys, Es, 
		range(start, end + 1),
		[neutrino], [kinematic_folder],
		charm_masses, 
		scales,
		parametrizations
	)

def differential(args, xs, ys, Es, neutrino, kinematic_folder):
	"""
	Generate an iterable of all argument combinations. Some arguments are read directly from command-line arguments.

	Args:
		args: Instance of argparse
		xs: List of Bjorken-x values
		ys: List of inelasticity (y) values
		Es: List of beam energy (E) values
		neutrino: If true, do neutrino scattering. Otherwise, do antineutrino scattering
		kinematic_folder: Name of the kinematic folder, used for generating the output directory structure
	"""
	charm_masses = [None]
	if args.masses and args.mode == Mode.masses:
		charm_masses = args.masses
	elif args.charm_mass is not None:
		charm_masses = [args.charm_mass]

	# By default, no scale variations
	scales = [(1.0, 1.0, 1.0)]
	# If the run mode is scales, then include all scale variations
	if args.mode == Mode.scales7p or args.mode == Mode.scales17p:
		scales = list(scale_variations.generate_standard(args))

	parametrizations = [None]
	if args.mode == Mode.decay_function_errors and args.decay_parametrizations:
		parametrization_list = np.loadtxt(args.decay_parametrizations, delimiter=",")
		parametrizations = [tuple(p) for p in parametrization_list]

	return full_differential(
		args,
		xs, ys, Es, 
		args.pdf_start, args.pdf_end, 
		neutrino, 
		kinematic_folder, 
		charm_masses, 
		scales,
		parametrizations[args.decay_start:(args.decay_end + 1)]
	)