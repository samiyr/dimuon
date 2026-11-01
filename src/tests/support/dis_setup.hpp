#ifndef DM_TESTS_SUPPORT_DIS_SETUP_H
#define DM_TESTS_SUPPORT_DIS_SETUP_H

#include <array>
#include <unordered_map>
#include <vector>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/cross_section_observable.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/scale_dependence.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/dis/channel_factory.hpp>
#include <dimuon/dis/instance.hpp>
#include <dimuon/dis/partonic_channel.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/pdf/is_pdf_interface.hpp>

namespace dm::testing {
	auto make_dis(
		const scattering_process::process_type type, const mass_scheme_t massless, const mass_scheme_t massive, const double charm_mass,
		const is_scalar_pdf_interface auto &pdf
	) {
		const scattering_process process{type, constants::particles::proton, constants::particles::neutrino};
		std::array<double, constants::total_flavors> masses{};
		masses[flavor::to_index(flavor::c)] = masses[flavor::to_index(flavor::cbar)] = charm_mass;
		const parton_configuration configuration{std::vector{flavor::u, flavor::d, flavor::s, flavor::c}, masses};
		std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> channels;
		for (const perturbative_order pto : {perturbative_order::lo, perturbative_order::nlo}) {
			channels[pto] = dis::channel_factory::create_all_channels(pto, process, configuration, channel_configuration{massless}, channel_configuration{massive});
		}
		integration_parameters parameters;
		parameters.cquad.max_relative_error = 1e-10;
		return dis::instance{
			ckm_matrix::global_fit(), pdf, {{perturbative_order::lo, parameters}, {perturbative_order::nlo, parameters}}, process,
			scale_dependence::multiplicative(1.0), scale_dependence::multiplicative(1.0),
			true, cross_section_observable{}, perturbative_order::nlo, 0.0, channels
		};
	}
}

#endif
