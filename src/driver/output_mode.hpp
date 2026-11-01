#ifndef DM_EXE_OUTPUT_MODE_H
#define DM_EXE_OUTPUT_MODE_H

/// Output mode of the results
enum class output_mode {
	/// Prints the results to stdout
	standard_output,
	/// Writes the results to a file and prints its path
	file
};

#endif