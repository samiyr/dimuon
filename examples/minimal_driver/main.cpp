#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <dimuon/common/cross_section_observable.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/decay/decay_grid.hpp>
#include <dimuon/pdf/extrapolator/freeze_extrapolator.hpp>
#include <dimuon/pdf/interface/lha_interface.hpp>
#include <dimuon/sidis/channel_factory.hpp>
#include <dimuon/sidis/instance.hpp>

int main() {
	using namespace dm;
	using namespace dm::constants::particles;

	// Define the scattering process: proton-neutrino scattering (proton due to incoherent scattering)
	const scattering_process process{scattering_process::process_type::neutrino, proton, neutrino};
	// Define the kinematical point (x = 0.1, y = 0.324, E = 88.29)
	const trf_kinematics kinematics = trf_kinematics::y_E_beam(0.1, 0.324, 88.29, proton.mass, 0.0);

	// Construct the PDF and FFs
	const lha_interface pdf{"EPPS21nlo_CT18Anlo_Fe56"};
	const std::vector ffs{
		lha_interface<freeze_extrapolator>{"kkks08_opal_d0___mas"},
		lha_interface<freeze_extrapolator>{"kkks08_opal_d+___mas"},
		lha_interface<freeze_extrapolator>{"bkk05_D3_d_s_nlo"},
		1.14 * lha_interface<freeze_extrapolator>{"bkk05_D3_lambda_c_nlo"}
	};

	// Define the decay muon energy cut
	const double muon_min_energy = 5.0;
	// Get the decay grid directory from an environment variable (the current directory if the variable is not set)
	const char *decay_grids_variable = std::getenv("DIMUON_DECAY_GRIDS_DIR");
	const std::filesystem::path decay_grids = decay_grids_variable != nullptr ? decay_grids_variable : "";
	// A helper lambda for constructing a decay instance for each intermediate resonance particle
	const auto decay = [&](const particle &hadron) {
		// The decay grid files follow a specific naming convention, which is defined by `decay_grid::filename`
		const std::string filename = decay_grid::filename(muon_min_energy, decay_parametrization::default_fit(), hadron, proton, massless_muon);
		const std::filesystem::path grid_path = decay_grids / filename;
		// A missing grid usually means that `DIMUON_DECAY_GRIDS_DIR` is not set or points to the wrong directory
		if (!std::filesystem::exists(grid_path)) {
			throw std::runtime_error(std::format("Decay grid {} not found. Did you set DIMUON_DECAY_GRIDS_DIR correctly?", grid_path.generic_string()));
		}
		// Construct the decay instance
		return decay_instance{hadron, proton, decay_grid{grid_path}, muon_min_energy};
	};
	// Define the list of decay instances
	const std::vector decays{decay(D0), decay(Dp), decay(Ds), decay(LambdaC)};
	// Define the allowed quark flavors
	const parton_configuration quarks{std::vector{parton{flavor::d}, parton{flavor::s}, parton{flavor::c, pdf.quark_mass(flavor::c)}}};
	// Construct all leading-order channels in the slow-rescaling scheme
	const auto channels = sidis::channel_factory::create_all_channels(
		perturbative_order::lo, process, quarks, std::nullopt, channel_configuration{mass_scheme_t::slow_rescaling}
	);
	// Construct the SIDIS instance
	const sidis::instance sidis{
		// Define the coupling (this is charged-current scattering, so use the CKM matrix elements)
		ckm_matrix::global_fit(),
		// Pass on the PDF, FFs, and decay instances
		pdf, ffs, decays,
		// Define the integration parameters (which are per perturbative order)
		{{perturbative_order::lo, integration_parameters{}}},
		// Pass the scattering process
		process,
		// Define the scale dependence functions (renormalization, factorization, fragmentation)
		// The `Q2` scale dependence is just Q^2
		scale_dependence::Q2, scale_dependence::Q2, scale_dependence::Q2,
		// Use a modified cross section prefactor to match the literature results
		true, 
		// Define the observable to be computed, in this case the neutrino cross section
		cross_section_observable{},
		// Define the highest perturbative order to be calculated
		perturbative_order::lo,
		// Define the minimum primary muon energy (set to zero, a non-zero value introduces a cut on `y`)
		0.0,
		// Define the partonic channels (which are per perturbative order)
		{{perturbative_order::lo, channels}},
		// Define the coefficient implementation, either `hw` or `gc` (which are per perturbative order)
		// The coefficient implementations are explained in more detail in `numerical_evaluation.md`
		{{perturbative_order::lo, coefficient_implementation_t::hw}},
		// Define the heavy-quark thresholds (which are set to zero here)
		heavy_quark_thresholds{}
	};

	// Calculate the cross section
	const sidis::result result = sidis.lepton_pair_cross_section_xy(kinematics);
	// Print the LO value
	std::cout << result.value.lo() << std::endl;
}
