import argparse
import json
from pathlib import Path

from pineappl_conversion_functions import dump_grid_directory_to_pineappl

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"-i", "--inputdir", type=Path, required=True, 
		help="Path to directory containing all grid files that should be converted"
	)
	parser.add_argument(
		"-o", "--outputdir", type=Path, required=True, 
		help="Path to directory where the converted grids will be written, in a subdirectory for each grid file"
	)
	parser.add_argument(
		"-r", "--resonances", type=json.loads, required=False, 
		default='{"D0": 1.0, "Dp": 1.0, "Ds": 1.0, "LambdaC": 1.14}', 
		help="Dictionary defining the resonance particles and their associated normalization constants. Names must match the ones in the input grids."
	)
	parser.add_argument(
		"-p", "--pattern", type=str, required=False, default="*.csv.gz", 
		help="Glob pattern for searching the inputdir"
	)

	args = parser.parse_args()

	dump_grid_directory_to_pineappl(args.inputdir, args.outputdir, args.resonances, args.pattern)