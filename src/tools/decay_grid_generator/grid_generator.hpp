#ifndef DM_TOOLS_DECAY_GRID_GENERATOR_GRID_GENERATOR_H
#define DM_TOOLS_DECAY_GRID_GENERATOR_GRID_GENERATOR_H

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <vector>
#include <string>
#include <format>
#include <fstream>

#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/decay/decay_grid.hpp>

#include <dimuon/common/particle.hpp>

#include <dimuon/utility/io.hpp>
#include <dimuon/utility/math.hpp>

#include <dimuon/integration/integrator_nd.hpp>
#include <dimuon/integration/integrator_cquad.hpp>

#include "decay_function.hpp"

struct grid_generator {
	template <typename T>
	static unsigned int number_of_digits(const T x) {
		return static_cast<unsigned int>(std::log10(static_cast<double>(x)) + 1.0);
	}

	static void generate_decay_grids(
		const std::filesystem::path grid_folder, const std::vector<double> &log_zyE_bins,
		const std::vector<double> &E_min_bins, const std::vector<dm::decay_parametrization> &parametrization_bins, 
		const std::vector<dm::particle> &resonance_bins, const dm::particle &target, const dm::particle &lepton
	) {
		const std::size_t member_count = E_min_bins.size() * parametrization_bins.size() * resonance_bins.size();

		#pragma omp parallel for collapse(3)
		for (std::size_t parametrization_index = 0; parametrization_index < parametrization_bins.size(); parametrization_index++) {
			for (std::size_t resonance_index = 0; resonance_index < resonance_bins.size(); resonance_index++) {
				for (std::size_t Emin_index = 0; Emin_index < E_min_bins.size(); Emin_index++) {
					const dm::particle &resonance = resonance_bins[resonance_index];
					const dm::decay_parametrization &parametrization = parametrization_bins[parametrization_index];
					const double E_min = E_min_bins[Emin_index];

					const std::string filename = dm::decay_grid::filename(E_min, parametrization, resonance, target, lepton);

					const std::filesystem::path grid_path = grid_folder / filename;

					dm::io::create_directory_tree(grid_path);
					std::ofstream file(grid_path);

					// The metadata, the grid points and the grid values, separated by "-----" (see dm::decay_grid, which reads the file)
					file << "version: " << dm::decay_grid::version << std::endl;
					file << "E_min: " << E_min << std::endl;
					file << "resonance: " << resonance << std::endl;
					file << "target: " << target << std::endl;
					file << "lepton: " << lepton << std::endl;
					file << "N: " << parametrization.N << std::endl;
					file << "alpha: " << parametrization.alpha << std::endl;
					file << "beta: " << parametrization.beta << std::endl;
					file << "gamma: " << parametrization.gamma << std::endl;
					file << "grid points: log(zyE)" << std::endl;
					file << "-----" << std::endl;

					for (const double value : log_zyE_bins) {
						file << std::format("{:.8f}", value) << std::endl;
					}
					file << "-----" << std::endl;

					const std::size_t zyE_count = log_zyE_bins.size();
					const std::size_t total_count = zyE_count;

					unsigned int calculated_values = 0;

					std::vector<double> grid(total_count);

					for (std::size_t i = 0; i < zyE_count; i++) {
						const double zyE = std::exp(log_zyE_bins[i]);

						const dm::integrator_nd<dm::integrator_cquad, 2> integrator(0.0, 1e-6);

						// Enforces the cut `rho E_h >= E_min` where `E_h = zyE`
						const double rho_min = std::min(E_min / zyE, 1.0);

						const double value = integrator.integrate(
							[&](const double rho, const double cos) {
								const double h0 = zyE;
								if (h0 < resonance.mass) { return 0.0; }

								const double mu = lepton.mass / h0;
								const double reduced_rho = sqrt(dm::math::power<2>(rho) - dm::math::power<2>(mu));

								const double a = dm::math::power<2>(h0) / dm::math::power<2>(resonance.mass);
								const double b = h0 * sqrt(dm::math::power<2>(h0) - dm::math::power<2>(resonance.mass)) / dm::math::power<2>(resonance.mass);

								return differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
							}, {rho_min, -1.0}, {1.0, 1.0}
						);

						grid[i] = value;

						calculated_values++;

						const int base_precision = 5;
						const int zyE_precision = base_precision - static_cast<int>(number_of_digits(static_cast<int>(zyE)));
						// The position of this grid in the loops, from 1 to member_count
						const std::size_t member_index = (parametrization_index * resonance_bins.size() + resonance_index) * E_min_bins.size() + Emin_index + 1;

						#pragma omp critical
						{
							std::cout << std::fixed << std::setprecision(base_precision);
							std::cout << "[grid] " << dm::io::leading_zeroes(calculated_values, number_of_digits(total_count));
							std::cout << " / " << total_count;
							std::cout << " [" << "grid " << dm::io::leading_zeroes(member_index, number_of_digits(member_count));
							std::cout << " / " << member_count << "]";
							std::cout << ": " << std::format("{:.8f}", value);
							std::cout << std::setprecision(zyE_precision) << " (zyE = " << zyE << ")";
							std::cout << "\r" << std::flush;
						}
					}

					for (std::size_t i = 0; i < zyE_count; i++) {
						file << std::format("{:.8f}", grid[i]) << std::endl;
					}

					file.close();
				}
			}
		}
	}
};

#endif