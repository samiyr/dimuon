#include <vector>								// std::vector
#include <string>								// std::string
#include <optional>								// std::optional
#include <iostream>								// std::cout, std::cerr, std::endl
#include <cstdlib>								// EXIT_SUCCESS, EXIT_FAILURE
#include <exception>							// std::exception

#include <argparser/arg_parser.hpp>				// arg_parser

#include <dimuon/common/constants.hpp>		// dm::constants::particles
#include <dimuon/utility/math.hpp>			// dm::math::chebyshev_space

#include "grid_generator.hpp"

std::vector<dm::decay_parametrization> read_input_file(const std::string input) {
	std::ifstream file(input);
	if (!file) {
		throw std::runtime_error("Fit set file '" + input + "' could not be opened.");
	}

	std::vector<std::vector<double>> parameter_set;

	std::size_t current_index = 0;

	std::string line;
	while (std::getline(file, line)) {
		current_index++;

		std::vector<double> parameters;

		std::istringstream stream(line);
		std::string value;

		while (std::getline(stream, value, ',')) {
			parameters.push_back(std::stod(value));
		}

		parameter_set.push_back(parameters);
	}

	std::vector<dm::decay_parametrization> parametrizations;

	for (const std::vector<double> &parameters : parameter_set) {
		parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
	}

	return parametrizations;
}

int run(const arg_parser &args) {
	using namespace dm;

	const double minE = args.get_required_double_argument("minE");

	// The grid starts by default at the minimum muon energy, below which the decay function vanishes
	const double bin_min = args.get_double_argument("bin_min", minE);
	const double bin_max = args.get_double_argument("bin_max", 300.0);
	const unsigned int bin_count = args.get_unsigned_int_argument("bin_count", 30);

	const double lepton_mass = args.get_double_argument("lepton_mass", 0.0);

	const std::optional<std::string> output_dir_opt = args.get_argument("output");
	const std::filesystem::path output_path = output_dir_opt ? std::filesystem::path(*output_dir_opt) : std::filesystem::current_path();
	const std::string output_dir = output_path.string() + "/";

	const std::string input_file = args.get_required_argument("input");

	// All arguments have been read
	args.check_arguments();

	std::cout << "Output directory: " << output_dir << std::endl;

	const std::vector<dm::decay_parametrization> parametrizations = read_input_file(input_file);
	// The grid points are Chebyshev points in log(zyE), the variable that dm::decay_grid interpolates in
	const std::vector<double> bins = math::chebyshev_space(std::log(bin_min), std::log(bin_max), bin_count);

	grid_generator::generate_decay_grids(
		output_dir + "decaygrids", bins, {minE}, parametrizations,
		{
			constants::particles::D0, constants::particles::Dp, constants::particles::Ds, constants::particles::LambdaC
		}, constants::particles::proton, particle(lepton_mass)
	);

	return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
	try {
		return run(arg_parser(argc, argv));
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}
