import numpy as np
import pandas as pd
from utils import launch
from utils import pool
import time
import datetime
import sys

if __name__ == '__main__':
	Wsign = sys.argv[1]
	output_file = sys.argv[2]

	run_args = [[z, z + 0.02, 25.0, 159.0, Wsign] for z in np.linspace(0.2, 0.84, 33)]
	run_args += [[z, z + 0.02, 159.0, 1007.0, Wsign] for z in np.linspace(0.2, 0.84, 33)]
	run_args += [[z, z + 0.02, 1007.0, 6400.0, Wsign] for z in np.linspace(0.2, 0.84, 33)]

	start = time.time()
	outputs = pool.run(launch.pion_production, run_args)

	# Create a dataframe from output list
	print("Creating dataframe")
	df = pd.DataFrame(outputs, columns=["z_min", "z_max", "Q2_min", "Q2_max", "LO", "NLO", "NNLO"])
	print("Dataframe shape:", df.shape)

	df.to_csv(output_file, index=False, header=True)

	delta_seconds = time.time() - start
	delta = datetime.timedelta(seconds=int(delta_seconds))
	print("Done, took", str(delta))