#include <cstddef>
#include <format>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include <dimuon/common/channel_configuration.hpp>
#include <dimuon/common/coefficient_implementation.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/cross_section_observable.hpp>
#include <dimuon/common/mass_scheme.hpp>
#include <dimuon/common/perturbative_order.hpp>
#include <dimuon/common/scale_dependence.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/dis/channel_factory.hpp>
#include <dimuon/dis/instance.hpp>
#include <dimuon/dis/partonic_channel.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/sidis/instance.hpp>
#include <dimuon/sidis/partonic_channel.hpp>
#include <dimuon/utility/io.hpp>

#include "support/sidis_setup.hpp"

using namespace dm;
using namespace dm::testing;

namespace {
	using enum perturbative_order;
	using enum coefficient_implementation_t;
	using enum scattering_process::process_type;

	const scattering_process neutrino_process{neutrino, constants::particles::proton, constants::particles::neutrino};

	/// A SIDIS configuration of neutrino scattering up to NNLO that is valid as given: that of `compute_sidis`, with the toy PDF and FF.
	/// The tests make one part of it invalid at a time.
	struct sidis_configuration {
		perturbative_order order = nnlo;
		std::unordered_map<perturbative_order, std::vector<sidis::partonic_channel>> channels = make_sidis_channels(sidis_run{});
		std::unordered_map<perturbative_order, coefficient_implementation_t> implementations{{lo, hw}, {nlo, hw}, {nnlo, gc}};
		std::unordered_map<perturbative_order, integration_parameters> integration{{lo, {}}, {nlo, {}}, {nnlo, {}}};
		std::size_t ff_count = 2;
		std::size_t decay_count = 2;

		/// Constructs the instance, which checks the configuration
		void construct() const {
			static io::null_stream progress;
			const decay_instance decay{test_resonances[0], constants::particles::proton, unit_decay_function{}, test_lepton_momentum_min};
			static_cast<void>(sidis::instance{
				ckm_matrix::global_fit(), toy_pdf_interface(), std::vector(ff_count, toy_ff_interface()), std::vector(decay_count, decay), integration, neutrino_process,
				scale_dependence::multiplicative(1.0), scale_dependence::multiplicative(1.0), scale_dependence::multiplicative(1.0),
				true, cross_section_observable{}, order, 0.0, channels, implementations, test_thresholds, progress
			});
		}
	};

	/// The DIS channels of the process \p type at LO and NLO, with the flavors and quark masses of the SIDIS tests (see `make_parton_configuration`)
	std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> make_dis_channels(const scattering_process::process_type type) {
		const scattering_process process{type, constants::particles::proton, constants::particles::neutrino};
		std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> channels;
		for (const perturbative_order pto : {lo, nlo}) {
			channels[pto] = dis::channel_factory::create_all_channels(
				pto, process, make_parton_configuration(1.3, 4.75), channel_configuration{mass_scheme_t::slow_rescaling}, channel_configuration{mass_scheme_t::sacot_subtracted}
			);
		}
		return channels;
	}

	/// A DIS configuration of neutrino scattering up to NLO that is valid as given, with the toy PDF. The tests make one part of it invalid at a time.
	struct dis_configuration {
		perturbative_order order = nlo;
		std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> channels = make_dis_channels(neutrino);
		std::unordered_map<perturbative_order, integration_parameters> integration{{lo, {}}, {nlo, {}}};

		/// Constructs the instance, which checks the configuration
		void construct() const {
			static_cast<void>(dis::instance{
				ckm_matrix::global_fit(), toy_pdf_interface(), integration, neutrino_process,
				scale_dependence::multiplicative(1.0), scale_dependence::multiplicative(1.0),
				true, cross_section_observable{}, order, 0.0, channels
			});
		}
	};
}

TEST(sidis_instance, valid_configuration) {
	EXPECT_NO_THROW(sidis_configuration{}.construct());

	// Up to NLO, nothing of NNLO is needed
	sidis_configuration up_to_nlo{.order = nlo};
	up_to_nlo.channels.erase(nnlo);
	up_to_nlo.implementations.erase(nnlo);
	up_to_nlo.integration.erase(nnlo);
	EXPECT_NO_THROW(up_to_nlo.construct());
}

TEST(sidis_instance, missing_orders) {
	// Each order up to the highest needs its channels, coefficient implementation, and integration parameters
	for (const perturbative_order pto : {lo, nlo, nnlo}) {
		sidis_configuration without_channels;
		without_channels.channels.erase(pto);
		EXPECT_THROW(without_channels.construct(), std::invalid_argument) << std::format("channels of {}", pto);

		sidis_configuration without_implementation;
		without_implementation.implementations.erase(pto);
		EXPECT_THROW(without_implementation.construct(), std::invalid_argument) << std::format("coefficient implementation of {}", pto);

		sidis_configuration without_integration;
		without_integration.integration.erase(pto);
		EXPECT_THROW(without_integration.construct(), std::invalid_argument) << std::format("integration parameters of {}", pto);
	}
}

TEST(sidis_instance, nnlo_hw) {
	// The NNLO coefficient functions have only the gc implementation, which is checked also when NNLO is not computed
	for (const perturbative_order highest_order : {nlo, nnlo}) {
		sidis_configuration configuration{.order = highest_order};
		configuration.implementations.at(nnlo) = hw;
		EXPECT_THROW(configuration.construct(), std::invalid_argument) << std::format("up to {}", highest_order);
	}
}

TEST(sidis_instance, sacot) {
	// At LO, the SACOT schemes are those of slow rescaling, also with gc (the massive channels of the configuration are in the SACOT scheme gm)
	sidis_configuration lo_gc;
	lo_gc.implementations.at(lo) = gc;
	EXPECT_NO_THROW(lo_gc.construct());

	// The mass terms of the SACOT schemes are only in the hw coefficient functions of NLO, which is checked also when NLO is not computed
	for (const perturbative_order highest_order : {lo, nlo}) {
		sidis_configuration nlo_gc{.order = highest_order};
		nlo_gc.implementations.at(nlo) = gc;
		EXPECT_THROW(nlo_gc.construct(), std::invalid_argument) << std::format("up to {}", highest_order);
	}

	// There are no SACOT coefficient functions at NNLO, which is checked also when NNLO is not computed
	for (const perturbative_order highest_order : {nlo, nnlo}) {
		const sidis_configuration nnlo_sacot{.order = highest_order, .channels = make_sidis_channels(sidis_run{.massive_scheme_nnlo = mass_scheme_t::sacot_subtracted})};
		EXPECT_THROW(nnlo_sacot.construct(), std::invalid_argument) << std::format("up to {}", highest_order);
	}
}

TEST(sidis_instance, invalid_channels) {
	// A channel listed with the channels of another order
	sidis_configuration misplaced;
	misplaced.channels.at(lo).push_back(misplaced.channels.at(nlo).front());
	EXPECT_THROW(misplaced.construct(), std::invalid_argument);

	// The channels of antineutrino scattering, which don't conserve the charge in neutrino scattering
	sidis_configuration antineutrino_channels;
	antineutrino_channels.channels = make_sidis_channels(sidis_run{.process = antineutrino});
	EXPECT_THROW(antineutrino_channels.construct(), std::invalid_argument);
}

TEST(sidis_instance, ffs_and_decays) {
	// Each decay channel needs an FF, and there must be at least one decay channel
	EXPECT_THROW((sidis_configuration{.ff_count = 1, .decay_count = 2}.construct()), std::invalid_argument);
	EXPECT_THROW((sidis_configuration{.ff_count = 2, .decay_count = 1}.construct()), std::invalid_argument);
	EXPECT_THROW((sidis_configuration{.ff_count = 0, .decay_count = 0}.construct()), std::invalid_argument);
}

TEST(dis_instance, valid_configuration) {
	EXPECT_NO_THROW(dis_configuration{}.construct());

	// At LO, nothing of NLO is needed
	dis_configuration lo_only{.order = lo};
	lo_only.channels.erase(nlo);
	lo_only.integration.erase(nlo);
	EXPECT_NO_THROW(lo_only.construct());
}

TEST(dis_instance, order_above_nlo) {
	// There are DIS coefficient functions only up to NLO, so NNLO is invalid even with its channels and integration parameters
	dis_configuration configuration{.order = nnlo};
	configuration.channels[nnlo] = {};
	configuration.integration[nnlo] = {};
	EXPECT_THROW(configuration.construct(), std::invalid_argument);
}

TEST(dis_instance, missing_orders) {
	// Each order up to the highest needs its channels and integration parameters
	for (const perturbative_order pto : {lo, nlo}) {
		dis_configuration without_channels;
		without_channels.channels.erase(pto);
		EXPECT_THROW(without_channels.construct(), std::invalid_argument) << std::format("channels of {}", pto);

		dis_configuration without_integration;
		without_integration.integration.erase(pto);
		EXPECT_THROW(without_integration.construct(), std::invalid_argument) << std::format("integration parameters of {}", pto);
	}
}

TEST(dis_instance, invalid_channels) {
	// A channel listed with the channels of another order
	dis_configuration misplaced;
	misplaced.channels.at(lo).push_back(misplaced.channels.at(nlo).front());
	EXPECT_THROW(misplaced.construct(), std::invalid_argument);

	// The channels of antineutrino scattering, which don't conserve the charge in neutrino scattering
	dis_configuration antineutrino_channels;
	antineutrino_channels.channels = make_dis_channels(antineutrino);
	EXPECT_THROW(antineutrino_channels.construct(), std::invalid_argument);
}
