#ifndef DM_TESTS_SUPPORT_SIDIS_SETUP_H
#define DM_TESTS_SUPPORT_SIDIS_SETUP_H

#include <array>
#include <cmath>
#include <flat_set>
#include <ranges>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <dimuon/common/coefficient_implementation.hpp>
#include <dimuon/common/constants.hpp>
#include <dimuon/common/cross_section_observable.hpp>
#include <dimuon/common/heavy_quark_thresholds.hpp>
#include <dimuon/common/parton_configuration.hpp>
#include <dimuon/common/scale_dependence.hpp>
#include <dimuon/common/scattering_process.hpp>
#include <dimuon/common/structure_function.hpp>
#include <dimuon/common/trf_kinematics.hpp>
#include <dimuon/coupling/ckm_matrix.hpp>
#include <dimuon/decay/decay_instance.hpp>
#include <dimuon/integration/integration_parameters.hpp>
#include <dimuon/integration/integrator_mc.hpp>
#include <dimuon/pdf/interface/functional_form_interface.hpp>
#include <dimuon/sidis/channel_factory.hpp>
#include <dimuon/sidis/evaluator.hpp>
#include <dimuon/sidis/instance.hpp>
#include <dimuon/utility/io.hpp>
#include <dimuon/utility/math.hpp>

namespace dm::testing {
	inline double toy_pdf(const flavor::type flavor, const double x) {
		const double i = static_cast<double>(flavor::to_index(flavor));
		return (1.0 + 0.1 * i) * std::pow(x, 0.3 + 0.07 * i) * std::pow(1.0 - x, 2.0 + 0.3 * i);
	}

	inline double toy_ff(const flavor::type flavor, const double z) {
		const double i = static_cast<double>(flavor::to_index(flavor));
		return (0.5 + 0.05 * i) * std::pow(z, 1.0 + 0.11 * i) * std::pow(1.0 - z, 1.5 + 0.2 * i);
	}

	inline auto toy_pdf_interface(const bool charge_conjugate = false) {
		return functional_form_interface{[charge_conjugate](const flavor::type f, const double x, const double) {
			return toy_pdf(charge_conjugate ? flavor::conjugate_flavor(f) : f, x);
		}};
	}

	inline auto toy_ff_interface(const bool charge_conjugate = false) {
		return functional_form_interface{[charge_conjugate](const flavor::type f, const double z, const double) {
			return toy_ff(charge_conjugate ? flavor::conjugate_flavor(f) : f, z);
		}};
	}

	inline const std::vector<flavor::type> sidis_flavors{flavor::u, flavor::d, flavor::s, flavor::c, flavor::b};

	inline parton_configuration make_parton_configuration(const double charm_mass, const double bottom_mass) {
		std::array<double, constants::total_flavors> masses{};
		for (const flavor::type f : {flavor::c, flavor::cbar}) { masses[flavor::to_index(f)] = charm_mass; }
		for (const flavor::type f : {flavor::b, flavor::bbar}) { masses[flavor::to_index(f)] = bottom_mass; }
		return parton_configuration{sidis_flavors, masses};
	}

	decltype(auto) with_coefficient_functions(const perturbative_order pto, const coefficient_implementation_t implementation, auto &&evaluate) {
		using namespace sidis::coefficients;
		const bool hw = implementation == coefficient_implementation_t::hw;
		switch (pto) {
			case perturbative_order::lo:
				if (hw) {
					return evaluate(
						f1_hw::lo::quark_to_quark, f1_hw::lo::quark_to_gluon, f1_hw::lo::gluon_to_quark, f1_hw::lo::gluon_to_gluon,
						f2_hw::lo::quark_to_quark, f2_hw::lo::quark_to_gluon, f2_hw::lo::gluon_to_quark, f2_hw::lo::gluon_to_gluon,
						f3_hw::lo::quark_to_quark, f3_hw::lo::quark_to_gluon, f3_hw::lo::gluon_to_quark, f3_hw::lo::gluon_to_gluon
					);
				}
				return evaluate(
					f1_gc::lo::quark_to_quark, f1_gc::lo::quark_to_gluon, f1_gc::lo::gluon_to_quark, f1_gc::lo::gluon_to_gluon,
					f2_gc::lo::quark_to_quark, f2_gc::lo::quark_to_gluon, f2_gc::lo::gluon_to_quark, f2_gc::lo::gluon_to_gluon,
					f3_gc::lo::quark_to_quark, f3_gc::lo::quark_to_gluon, f3_gc::lo::gluon_to_quark, f3_gc::lo::gluon_to_gluon
				);
			case perturbative_order::nlo:
				if (hw) {
					return evaluate(
						f1_hw::nlo::quark_to_quark, f1_hw::nlo::quark_to_gluon, f1_hw::nlo::gluon_to_quark, f1_hw::nlo::gluon_to_gluon,
						f2_hw::nlo::quark_to_quark, f2_hw::nlo::quark_to_gluon, f2_hw::nlo::gluon_to_quark, f2_hw::nlo::gluon_to_gluon,
						f3_hw::nlo::quark_to_quark, f3_hw::nlo::quark_to_gluon, f3_hw::nlo::gluon_to_quark, f3_hw::nlo::gluon_to_gluon
					);
				}
				return evaluate(
					f1_gc::nlo::quark_to_quark, f1_gc::nlo::quark_to_gluon, f1_gc::nlo::gluon_to_quark, f1_gc::nlo::gluon_to_gluon,
					f2_gc::nlo::quark_to_quark, f2_gc::nlo::quark_to_gluon, f2_gc::nlo::gluon_to_quark, f2_gc::nlo::gluon_to_gluon,
					f3_gc::nlo::quark_to_quark, f3_gc::nlo::quark_to_gluon, f3_gc::nlo::gluon_to_quark, f3_gc::nlo::gluon_to_gluon
				);
			case perturbative_order::nnlo:
				if (hw) { throw std::invalid_argument("The NNLO coefficient functions have only the gc implementation"); }
				return evaluate(
					f1_gc::nnlo::quark_to_quark, f1_gc::nnlo::quark_to_gluon, f1_gc::nnlo::gluon_to_quark, f1_gc::nnlo::gluon_to_gluon,
					f2_gc::nnlo::quark_to_quark, f2_gc::nnlo::quark_to_gluon, f2_gc::nnlo::gluon_to_quark, f2_gc::nnlo::gluon_to_gluon,
					f3_gc::nnlo::quark_to_quark, f3_gc::nnlo::quark_to_gluon, f3_gc::nnlo::gluon_to_quark, f3_gc::nnlo::gluon_to_gluon
				);
		}
		std::unreachable();
	}

	struct structure_function_point {
		double x;
		double z;
		double Q2;
		double renormalization = 1.0;
		double factorization = 1.0;
		double fragmentation = 1.0;
	};

	inline constexpr heavy_quark_thresholds test_thresholds{1.3, 4.75, 172.5};

	struct unit_decay_function {
		double operator()(const double, const double, const double, const particle &) const noexcept { return 1.0; }
	};

	template <structure_function_t sf, is_scalar_pdf_interface pdf_t, is_scalar_pdf_interface ff_t, is_coupling coupling_t>
	class fixed_z_integrand {
	public:
		explicit fixed_z_integrand(
			const std::vector<sidis::partonic_channel> &channels_, const perturbative_order pto_, const coefficient_implementation_t implementation_,
			const pdf_t &pdf, const ff_t &ff, const coupling_t &coupling, const structure_function_point &point_
		) :
			pto{pto_},
			implementation{implementation_},
			point{point_},
			params{
				pdf, std::vector{ff}, coupling, process, trf_kinematics::Q2_s(point.x, point.Q2, 1000.0, process.target.mass, process.projectile.mass),
				pto, channels_, implementation, test_thresholds, point.renormalization * point.Q2, point.factorization * point.Q2, point.fragmentation * point.Q2
			},
			evaluator{channels_, channels_ | std::views::transform([](const sidis::partonic_channel &channel) { return sidis::physical_channel{channel}; }) | std::ranges::to<std::flat_set>()} {}

		double operator()(const double xi, const double zeta) const {
			double total = 0.0;
			with_coefficient_functions(pto, implementation, [&](auto &&...coefficient_functions) {
				evaluator.evaluate_observable(
					xi, zeta, point.z, params, decays, coefficient_functions..., structure_function::observable<sf>{},
					[&](const std::size_t, const std::size_t, const std::size_t, const double value) { total += value; }
				);
			});
			return total;
		}

	private:
		inline static const scattering_process process{scattering_process::process_type::neutrino, constants::particles::proton, constants::particles::neutrino};
		perturbative_order pto;
		coefficient_implementation_t implementation;
		structure_function_point point;
		sidis::parameters<pdf_t, ff_t, coupling_t> params;
		decay_parameters<unit_decay_function> decays{{decay_instance{particle{0.0, 1.0}, constants::particles::proton, unit_decay_function{}, 0.0}}};
		sidis::evaluator evaluator;
	};

	template <structure_function_t sf>
	double fixed_z_structure_function(
		const std::vector<sidis::partonic_channel> &channels, const perturbative_order pto, const coefficient_implementation_t implementation,
		const is_scalar_pdf_interface auto &pdf, const is_scalar_pdf_interface auto &ff, const is_coupling auto &coupling,
		const structure_function_point &point, const integration_parameters &parameters
	) {
		const fixed_z_integrand<sf, std::remove_cvref_t<decltype(pdf)>, std::remove_cvref_t<decltype(ff)>, std::remove_cvref_t<decltype(coupling)>> integrand{
			channels, pto, implementation, pdf, ff, coupling, point
		};
		return integrator_mc{parameters}.integrate([&](const std::span<const double> convolution, const std::span<double> components) {
			components[0] += integrand(convolution[0], convolution[1]);
		}, {point.x, 0.0}, {1.0, 1.0}, 1).value();
	}

	struct sidis_run {
		scattering_process::process_type process = scattering_process::process_type::neutrino;
		perturbative_order order = perturbative_order::nnlo;
		mass_scheme_t massless_scheme = mass_scheme_t::slow_rescaling;
		mass_scheme_t massive_scheme = mass_scheme_t::sacot_subtracted;
		mass_scheme_t massive_scheme_nnlo = mass_scheme_t::slow_rescaling;
		coefficient_implementation_t lo_implementation = coefficient_implementation_t::hw;
		coefficient_implementation_t nlo_implementation = coefficient_implementation_t::hw;
		double charm_mass = 1.3;
		double bottom_mass = 4.75;
		double x = 0.1;
		double y = 0.5;
		double E = 100.0;
		double renormalization = 1.0;
		double factorization = 1.0;
		double fragmentation = 1.0;
		integration_method method = integration_method::plain_qmc;
		long long points = 4096;
	};

	inline std::unordered_map<perturbative_order, std::vector<sidis::partonic_channel>> make_sidis_channels(const sidis_run &run) {
		const scattering_process process{run.process, constants::particles::proton, constants::particles::neutrino};
		const parton_configuration configuration = make_parton_configuration(run.charm_mass, run.bottom_mass);
		std::unordered_map<perturbative_order, std::vector<sidis::partonic_channel>> channels;
		for (const perturbative_order pto : {perturbative_order::lo, perturbative_order::nlo, perturbative_order::nnlo}) {
			const mass_scheme_t massive = pto == perturbative_order::nnlo ? run.massive_scheme_nnlo : run.massive_scheme;
			channels[pto] = sidis::channel_factory::create_all_channels(
				pto, process, configuration, channel_configuration{run.massless_scheme}, channel_configuration{massive}
			);
		}
		return channels;
	}

	inline double test_decay_function(const double z) { return 1.0 + z; }
	inline constexpr std::array test_resonances{particle{1.87, 1.0}, particle{2.29, 1.0}};
	inline constexpr double test_lepton_momentum_min = 2.0;

	sidis::result compute_sidis(const sidis_run &run, const is_pdf_interface auto &pdf, const is_pdf_interface auto &ff, const is_observable auto &observable) {
		const scattering_process process{run.process, constants::particles::proton, constants::particles::neutrino};
		const trf_kinematics kinematics = trf_kinematics::y_E_beam(run.x, run.y, run.E, process.target.mass, process.projectile.mass);

		integration_parameters parameters;
		parameters.method = run.method;
		parameters.cuba.maximum_evaluations = run.points;

		parameters.cuba.vegas.starting_iterations = parameters.cuba.suave.new_subdivision_evaluations = static_cast<int>(run.points / 10);
		const std::unordered_map<perturbative_order, integration_parameters> integration{
			{perturbative_order::lo, parameters}, {perturbative_order::nlo, parameters}, {perturbative_order::nnlo, parameters}
		};

		const auto decay_function = [](const double, const double z, const double, const particle &) { return test_decay_function(z); };
		const std::vector decays{
			decay_instance{test_resonances[0], constants::particles::proton, decay_function, test_lepton_momentum_min},
			decay_instance{test_resonances[1], constants::particles::proton, decay_function, test_lepton_momentum_min}
		};

		static io::null_stream progress;
		const sidis::instance instance{
			ckm_matrix::global_fit(), pdf, std::vector(decays.size(), ff), decays, integration, process,
			scale_dependence::multiplicative(run.renormalization), scale_dependence::multiplicative(run.factorization), scale_dependence::multiplicative(run.fragmentation),
			true, observable, run.order, 0.0, make_sidis_channels(run),
			std::unordered_map<perturbative_order, coefficient_implementation_t>{
				{perturbative_order::lo, run.lo_implementation}, {perturbative_order::nlo, run.nlo_implementation}, {perturbative_order::nnlo, coefficient_implementation_t::gc}
			},
			test_thresholds, progress
		};
		return instance.lepton_pair_cross_section_xy(kinematics);
	}

	sidis::result compute_sidis(const sidis_run &run, const is_pdf_interface auto &pdf, const is_pdf_interface auto &ff) {
		return compute_sidis(run, pdf, ff, cross_section_observable{});
	}
}

#endif
