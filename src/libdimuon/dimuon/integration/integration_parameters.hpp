#ifndef DM_INTEGRATION_PARAMETERS_H
#define DM_INTEGRATION_PARAMETERS_H

#include <array>
#include <cstddef>
#include <format>
#include <ostream>
#include <stdexcept>

#include "dimuon/utility/names.hpp"

namespace dm {

/// Integration method
enum class integration_method {
	/// Vegas (see `integrator_mc::integrate`)
	cuba_vegas,
	/// Suave (see `integrator_mc::integrate`)
	cuba_suave,
	/// Plain quasi-Monte Carlo with Sobol points (see `integrator_mc::integrate`)
	plain_qmc
};

/// Mapping of integration methods to their corresponding (argument) names, see `names`
inline constexpr auto integration_method_names = std::to_array<names::entry<integration_method>>({
	{integration_method::cuba_vegas, "vegas"},
	{integration_method::cuba_suave, "suave"},
	{integration_method::plain_qmc, "qmc"}
});

}

/// Formats an integration method as its name (see `dm::integration_method_names`)
template <> struct std::formatter<dm::integration_method> : dm::names::formatter<dm::integration_method_names> {};

namespace dm {

struct integration_parameters {
	struct cuba_parameters {
		struct vegas_parameters {
			int starting_iterations = 100'000;
			int continuing_iterations = 0;
			int batch_size = 1000;
			int grid_number = 0;

			friend std::ostream &operator <<(std::ostream &os, const vegas_parameters &v) {
				std::println(os, "#cuba.vegas.starting_iterations = {}", v.starting_iterations);
				std::println(os, "#cuba.vegas.continuing_iterations = {}", v.continuing_iterations);
				std::println(os, "#cuba.vegas.batch_size = {}", v.batch_size);
				std::println(os, "#cuba.vegas.grid_number = {}", v.grid_number);
				return os;
			}
		};

		struct suave_parameters {
			int new_subdivision_evaluations = 100'000;
			int minimum_samples = 25;
			double flatness = 50.0;

			friend std::ostream &operator <<(std::ostream &os, const suave_parameters &s) {
				std::println(os, "#cuba.suave.new_subdivision_evaluations = {}", s.new_subdivision_evaluations);
				std::println(os, "#cuba.suave.minimum_samples = {}", s.minimum_samples);
				std::println(os, "#cuba.suave.flatness = {}", s.flatness);
				return os;
			}
		};

		vegas_parameters vegas;
		suave_parameters suave;

		int seed = 0;
		long long int minimum_evaluations = 0;
		long long int maximum_evaluations = 1'000'000;

		double maximum_relative_error = 1e-3;
		double maximum_absolute_error = 1e-12;

		bool use_all_samples = true;
		bool apply_smoothing = true;
		bool delete_state_file_after_integration = true;
		bool reset_state = false;
		int random_number_generator_level = 0;

		bool parallelize = false;

		int verbose = 0;

		/// Constructs an integer representing the flag values in the format expected by Cuba
		int flags() const {
			int flag_value = 0;

			// Bits 0 and 1
			if (verbose < 0 || verbose > 3) {
				throw std::runtime_error(std::format("Cuba verbose flag must be an integer between 0 and 3, was given {} instead", verbose));
			}
			flag_value += verbose;
			// Bit 2
			flag_value += 4 * (!use_all_samples);
			// Bit 3
			flag_value += 8 * (!apply_smoothing);
			// Bit 4
			flag_value += 16 * (!delete_state_file_after_integration);
			// Bit 5
			flag_value += 32 * reset_state;
			// Bits 8-31
			flag_value += 256 * random_number_generator_level;

			return flag_value;
		}
	};

	struct cquad_parameters {
		std::size_t interval_count = 1'000;
		double max_absolute_error = 0.0;
		double max_relative_error = 1e-3;

		friend std::ostream &operator <<(std::ostream &os, const cquad_parameters &c) {
			std::println(os, "#cquad.interval_count = {}", c.interval_count);
			std::println(os, "#cquad.max_absolute_error = {}", c.max_absolute_error);
			std::println(os, "#cquad.max_relative_error = {}", c.max_relative_error);
			return os;
		}
	};

	integration_method method = integration_method::cuba_suave;
	cuba_parameters cuba;
	cquad_parameters cquad;

	friend std::ostream &operator <<(std::ostream &os, const integration_parameters &p) {
		std::println(os, "#method = {}", p.method);
		std::println(os, "#cuba.seed = {}", p.cuba.seed);
		std::println(os, "#cuba.minimum_evaluations = {}", p.cuba.minimum_evaluations);
		std::println(os, "#cuba.maximum_evaluations = {}", p.cuba.maximum_evaluations);

		std::println(os, "#cuba.maximum_relative_error = {}", p.cuba.maximum_relative_error);
		std::println(os, "#cuba.maximum_absolute_error = {}", p.cuba.maximum_absolute_error);

		// The booleans are printed as 1 or 0
		std::println(os, "#cuba.use_all_samples = {:d}", p.cuba.use_all_samples);
		std::println(os, "#cuba.apply_smoothing = {:d}", p.cuba.apply_smoothing);
		std::println(os, "#cuba.delete_state_file_after_integration = {:d}", p.cuba.delete_state_file_after_integration);
		std::println(os, "#cuba.reset_state = {:d}", p.cuba.reset_state);
		std::println(os, "#cuba.random_number_generator_level = {}", p.cuba.random_number_generator_level);

		os << p.cuba.vegas;
		os << p.cuba.suave;
		os << p.cquad;

		return os;
	}
};

}

#endif
