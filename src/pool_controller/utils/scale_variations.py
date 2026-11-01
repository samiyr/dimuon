import itertools
from .arg_modes import Mode

def generate_standard(args):
	"""
	Generate scale variations based on the command-line arguments.
	For SIDIS, all either three scale variations are included (renormalization, factorization, fragmentation) [scales17p]
	or the fragmentation scale is fixed to 1 [scales7p]. The extreme variations are removed. 
	For DIS, only renormalization and factorization are varied as there
	is no fragmentation scale. For technical simplicity, DIS scale variations still include the fragmentation scale,
	but it's fixed to 1. This keeps the tuple structure the same.

	Args:
		args: Arguments parsed by arg_parser
	"""
	scale_variations = [0.25, 1.0, 4.0]
	min_scale = min(scale_variations)
	max_scale = max(scale_variations)

	repeats = {
		Mode.scales7p: 2,
		Mode.scales17p: 3
	}
	repeat = repeats[args.mode]

	# Create all combinations of scale variations
	all_combinations = itertools.product(scale_variations, repeat=repeat)
	# If only renormalization and factorization scales are needed, fix fragmentation to 1
	if repeat == 2:
		all_combinations = [(renormalization, factorization, 1.0) for (renormalization, factorization) in all_combinations]

	def scale_filter_17p(arg):
		(renormalization, factorization, fragmentation) = arg
		factorization_condition = min_scale <= factorization / renormalization and max_scale >= factorization / renormalization
		fragmentation_condition = min_scale <= fragmentation / renormalization and max_scale >= fragmentation / renormalization

		return not (factorization_condition and fragmentation_condition)

	def scale_filter_7p(arg):
		(renormalization, factorization, _) = arg
		factorization_condition = min_scale <= factorization / renormalization and max_scale >= factorization / renormalization

		return not factorization_condition

	scale_filters = {
		Mode.scales7p: scale_filter_7p,
		Mode.scales17p: scale_filter_17p
	}

	return itertools.filterfalse(scale_filters[args.mode], all_combinations)