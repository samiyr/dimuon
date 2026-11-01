#ifndef DM_EXE_GRIDS_H
#define DM_EXE_GRIDS_H

#include <filesystem>
#include <array>
#include <utility>
#include <vector>
#include <cstdlib>
#include <format>
#include <stdexcept>
#include <string>
#include <system_error>

#include <dimuon/common/particle.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/decay/decay_parametrization.hpp>
#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/decay/decay_grid.hpp>
#include <dimuon/pdf/interface/lha_interface.hpp>
#include <dimuon/pdf/extrapolator/freeze_extrapolator.hpp>

namespace grids {
	/// @brief Returns the directory of the decay grids. It is the first one that exists of
	/// 1. the directory in the environment variable `DIMUON_DECAY_GRIDS_DIR`,
	/// 2. `share/dimuon/decaygrids` next to the `bin` directory of this program, which is where `cmake --install` puts the grids, and
	/// 3. `decaygrids` in the current directory.
	/// @throws std::runtime_error if none of them exists
	inline std::filesystem::path decay_grids_directory() {
		std::vector<std::filesystem::path> candidates;

		const char *directory = std::getenv("DIMUON_DECAY_GRIDS_DIR");
		if (directory && *directory) {
			candidates.emplace_back(directory);
		}

		// The file of this program is found through /proc, which only exists on Linux. Elsewhere this candidate is skipped.
		std::error_code error;
		const std::filesystem::path program = std::filesystem::read_symlink("/proc/self/exe", error);
		if (!error) {
			candidates.push_back(program.parent_path().parent_path() / "share" / "dimuon" / "decaygrids");
		}

		candidates.push_back(std::filesystem::current_path() / "decaygrids");

		std::string tried;
		for (const std::filesystem::path &candidate : candidates) {
			if (std::filesystem::is_directory(candidate)) {
				return candidate;
			}
			tried += std::format(" '{}'", candidate.string());
		}
		throw std::runtime_error(std::format(
			"No directory of decay grids exists (set the environment variable DIMUON_DECAY_GRIDS_DIR to one). Tried:{}", tried
		));
	}

	/// @brief Returns a function that constructs the decay of a single resonance, as a list of decays with one element
	/// @param[in] resonance The resonance
	inline auto create_resonance_decays(const dm::particle resonance) {
		using namespace dm;

		const auto resonance_decays = [=](
			const double decay_muon_min_energy, const particle &decay_lepton,
			const decay_parametrization &parametrization, const scattering_process &
		) {
			const std::filesystem::path decay_grid_folder = decay_grids_directory();
			const std::string decay_grid_filename = decay_grid::filename(
				decay_muon_min_energy, parametrization, resonance, constants::particles::proton, decay_lepton
			);
			const decay_grid decay{decay_grid_folder / decay_grid_filename};

			return std::vector{decay_instance{resonance, constants::particles::proton, decay, decay_muon_min_energy}};
		};

		return resonance_decays;
	}

	/// @brief Returns \p fragmentation_flavors without the fragmenting flavor that is masked in \p process, whose FFs are taken as zero:
	/// the charm antiquark in neutrino scattering, and the charm quark in antineutrino scattering.
	/// The channels are then restricted to the remaining flavors (see `channels::filter_channels`), so those of the masked flavor are not computed.
	inline dm::flavor::flavor_vector mask_fragmentation_flavors(dm::flavor::flavor_vector fragmentation_flavors, const dm::scattering_process &process) {
		using namespace dm;

		switch(process.type) {
			case scattering_process::process_type::neutrino:
				std::erase(fragmentation_flavors, flavor::cbar);
				break;
			case scattering_process::process_type::antineutrino:
				std::erase(fragmentation_flavors, flavor::c);
				break;
			default: break;
		}
		return fragmentation_flavors;
	}

	/// @brief Returns a function that constructs the FFs and the decays of the D0, D+, Ds and LambdaC resonances, in the same order.
	/// The FFs are those of the KKKS08 (D0, D+) and BKK05 (Ds, LambdaC) sets, with the members \p ff_members.
	inline auto create_grid_fragmentation(const std::array<std::size_t, 4> ff_members) {
		using namespace dm;

		const auto construct_grid_fragmentation_configuration = [=](
			const double decay_muon_min_energy, const decay_parametrization &parametrization,
			const particle &target, const particle &decay_lepton
		) {
			const std::filesystem::path decay_grid_folder = decay_grids_directory();
			const std::string D0_decay_grid = decay_grid::filename(
				decay_muon_min_energy, parametrization, constants::particles::D0, target, decay_lepton
			);
			const std::string Dp_decay_grid = decay_grid::filename(
				decay_muon_min_energy, parametrization, constants::particles::Dp, target, decay_lepton
			);
			const std::string Ds_decay_grid = decay_grid::filename(
				decay_muon_min_energy, parametrization, constants::particles::Ds, target, decay_lepton
			);
			const std::string LambdaC_decay_grid = decay_grid::filename(
				decay_muon_min_energy, parametrization, constants::particles::LambdaC, target, decay_lepton
			);

			const decay_grid D0_decay(decay_grid_folder / D0_decay_grid);
			const decay_grid Dp_decay(decay_grid_folder / Dp_decay_grid);
			const decay_grid Ds_decay(decay_grid_folder / Ds_decay_grid);
			const decay_grid LambdaC_decay(decay_grid_folder / LambdaC_decay_grid);

			// Halve the light-quark flavors when using KKKS08 / BKK05 FFs
			const std::array<double, constants::total_flavors> flavor_multipliers{1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

			// The FFs and the decays of the resonances, in the same order
			return std::pair{
				std::vector{
					lha_interface<freeze_extrapolator>("kkks08_opal_d0___mas", ff_members[0], flavor_multipliers),
					lha_interface<freeze_extrapolator>("kkks08_opal_d+___mas", ff_members[1], flavor_multipliers),
					lha_interface<freeze_extrapolator>("bkk05_D3_d_s_nlo", ff_members[2], flavor_multipliers),
					1.14 * lha_interface<freeze_extrapolator>("bkk05_D3_lambda_c_nlo", ff_members[3], flavor_multipliers)
				},
				std::vector{
					decay_instance(constants::particles::D0, target, D0_decay, decay_muon_min_energy),
					decay_instance(constants::particles::Dp, target, Dp_decay, decay_muon_min_energy),
					decay_instance(constants::particles::Ds, target, Ds_decay, decay_muon_min_energy),
					decay_instance(constants::particles::LambdaC, target, LambdaC_decay, decay_muon_min_energy)
				}
			};
		};

		const auto grid_fragmentation = [=](
			const double decay_muon_min_energy, const particle &decay_lepton, 
			const decay_parametrization &parametrization, const scattering_process &
		) {
			return construct_grid_fragmentation_configuration(
				decay_muon_min_energy, 
				parametrization, 
				constants::particles::proton, 
				decay_lepton
			);
		};

		return grid_fragmentation;
	}
}

#endif