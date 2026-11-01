def run(function, run_args, mpi:bool = False) -> list:
	"""
	Calls a function on each element of run_args in a thread pool, and returns the results in the same order.

	Args:
		function: Function called with one element of run_args, e.g. launch.dimuon_cross_section. It has to be defined at the top level of a module, since the processes of the pool look it up by name.
		run_args: Arguments of the computations, one element per call
		mpi: Whether to use an MPI pool (mpi4py.futures) instead of a local thread pool (multiprocessing)
	"""
	if mpi:
		from mpi4py.futures import MPIPoolExecutor
		# Initialize MPI pool
		print("Initializing MPI pool")
		executor = MPIPoolExecutor()
		# Iterate over all run arguments
		print("Starting iteration")
		outputs = executor.map(function, run_args)
		# Clean the pool after tasks are completed
		executor.shutdown()
		print("Pool shut down")
	else:
		from multiprocessing import Pool
		# Initialize pool
		print("Initializing pool")
		executor = Pool()
		# Iterate over all run arguments
		print("Starting iteration")
		outputs = executor.map(function, run_args)
		# Clean the pool after tasks are completed
		executor.close()
		executor.terminate()
		print("Pool shut down")

	return list(outputs)
