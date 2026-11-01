#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <optional>
#include <ranges>
#include <iostream>
#include <unordered_map>
#include <iomanip>
#include <cstdlib>
#include <concepts>
#include <filesystem>
#include <fstream>

#include <argparser/arg_parser.hpp>
#include "channels.hpp"
#include "maps.hpp"
#include "grids.hpp"
#include "run_mode.hpp"
#include "output_mode.hpp"

#include <dimuon/pdf/interface/vector_kernel_interface.hpp>
#include <dimuon/sidis/partonic_channel.hpp>

#include <dimuon/dis/instance.hpp>
#include <dimuon/sidis/instance.hpp>

#include <dimuon/coupling/ckm_matrix.hpp>

#include <dimuon/interpolation/linear_interpolation_basis.hpp>

#include <dimuon/utility/io.hpp>

static bool print_debug = true;

template <typename channel_t>
void print_channels(
	const std::unordered_map<dm::perturbative_order, std::vector<channel_t>> &channels, 
	const std::unordered_map<dm::perturbative_order, dm::coefficient_implementation_t> &coefficient_implementation
) {
	if (!print_debug) { return; }
	const auto print_order = [&](const dm::perturbative_order order, const std::string label) {
		if (channels.at(order).empty()) { return; }
		std::cout << label << "\t[" << dm::names::name_of(dm::coefficient_implementation_names, coefficient_implementation.at(order)) << "]" << std::endl;
		for (const channel_t &channel : channels.at(order)) {
			std::cout << channel << std::endl;
		}
		std::cout << std::endl;
	};

	if (channels.contains(dm::perturbative_order::lo)) { print_order(dm::perturbative_order::lo, "LO channels"); }
	if (channels.contains(dm::perturbative_order::nlo)) { print_order(dm::perturbative_order::nlo, "NLO channels"); }
	if (channels.contains(dm::perturbative_order::nnlo)) { print_order(dm::perturbative_order::nnlo, "NNLO channels"); }
}

/// Throws if there are no channels at any perturbative order
template <typename channel_t>
void check_channels(const std::unordered_map<dm::perturbative_order, std::vector<channel_t>> &channels) {
	const auto empty = [](const std::vector<channel_t> &order_channels) { return order_channels.empty(); };
	if (std::ranges::all_of(channels | std::views::values, empty)) {
		throw std::runtime_error("No channels: none are given (e.g. --nlo-all-m, --nlo-ba-m), or none have the flavors of --initial-flavors and --fragmentation-flavors (the masked fragmenting flavor is left out, see --disable-fragmentation-channel-mask)");
	}
}

/// The highest perturbative order that has channels, or LO if no order has them (which `check_channels` rules out)
template <typename channel_t>
dm::perturbative_order highest_order(const std::unordered_map<dm::perturbative_order, std::vector<channel_t>> &channels) {
	dm::perturbative_order order = dm::perturbative_order::lo;
	for (const auto &[channels_order, order_channels] : channels) {
		if (!order_channels.empty()) { order = std::max(order, channels_order); }
	}
	return order;
}

/// Writes the results with `write(stream)`: to stdout in the stdout output mode, and in the file output mode to the file \p output_file,
/// or if it is not given, to a new file with a random name in $TMPDIR. The path of the file is then printed.
template <std::invocable<std::ostream &> write_t>
void write_output(const output_mode mode, const std::optional<std::string> &output_file, const write_t &write) {
	if (mode == output_mode::standard_output) {
		write(std::cout);
		return;
	}

	const std::filesystem::path path = output_file ? std::filesystem::path(*output_file) : dm::io::create_temporary_file(".csv");
	std::ofstream file(path);
	if (!file) {
		throw std::runtime_error(std::format("Could not open the output file '{}'", path.string()));
	}
	file << std::setprecision(12);
	write(file);
	file.close();
	if (!file) {
		throw std::runtime_error(std::format("Could not write the output file '{}'", path.string()));
	}

	std::cout << "Output file written to:" << std::endl << path.string() << std::endl;
}

/// Computes the observable that the arguments \p args specify
int run(const arg_parser &args) {
	using namespace dm;

	lha_interface<>::disable_verbosity();

	std::cout << std::setprecision(12);

	print_debug = !args.get_flag("silent");

	// Read the observable and the run mode
	const std::string observable_name = args.get_required_argument("observable");
	const run_mode mode = args.get_required_choice_argument("run-mode", maps::run_modes);
	if (mode == run_mode::kernel && observable_name != "sidis.differential") {
		throw std::runtime_error("The kernel mode is only available for the observable sidis.differential");
	}

	// Read the output mode: the results are printed to stdout (default), or written to a file whose path is printed (see `write_output`).
	// The file is --output-file, or a new file with a random name in $TMPDIR.
	const output_mode output = args.get_choice_argument("output-mode", maps::output_modes, output_mode::standard_output);
	const std::optional<std::string> output_file = output == output_mode::file ? args.get_argument("output-file") : std::nullopt;

	// Read integration parameters
	const integration_parameters default_integration_parameters;
	integration_parameters common_integration_parameters;
	common_integration_parameters.cuba.minimum_evaluations = args.get_long_long_int_argument(
		"cuba-min-evals", default_integration_parameters.cuba.minimum_evaluations
	);
	common_integration_parameters.cuba.maximum_evaluations = args.get_long_long_int_argument(
		"cuba-max-evals", default_integration_parameters.cuba.maximum_evaluations
	);
	common_integration_parameters.cuba.maximum_absolute_error = args.get_double_argument(
		"cuba-max-abs-error", default_integration_parameters.cuba.maximum_absolute_error
	);
	common_integration_parameters.cuba.maximum_relative_error = args.get_double_argument(
		"cuba-max-rel-error", default_integration_parameters.cuba.maximum_relative_error
	);
	// The iteration sizes of Suave and Vegas default to a tenth of the maximum evaluations, i.e. 10 iterations
	const int default_iteration_evaluations = static_cast<int>(common_integration_parameters.cuba.maximum_evaluations / 10);
	common_integration_parameters.cuba.suave.new_subdivision_evaluations = args.get_int_argument(
		"cuba-suave-subd-evals", default_iteration_evaluations
	);
	common_integration_parameters.cuba.suave.minimum_samples = args.get_int_argument(
		"cuba-suave-min-samples", default_integration_parameters.cuba.suave.minimum_samples
	);
	common_integration_parameters.cuba.suave.flatness = args.get_double_argument(
		"cuba-suave-flatness", default_integration_parameters.cuba.suave.flatness
	);
	// The integration method: Suave (default), Vegas, or plain quasi-Monte Carlo in the direct mode.
	// The kernel mode integrates with plain quasi-Monte Carlo.
	common_integration_parameters.method = mode == run_mode::direct
		? args.get_choice_argument("integration-method", maps::integration_methods, integration_method::cuba_suave)
		: integration_method::plain_qmc;
	common_integration_parameters.cuba.vegas.starting_iterations = args.get_int_argument(
		"cuba-vegas-start-evals", default_iteration_evaluations
	);
	common_integration_parameters.cuba.vegas.continuing_iterations = args.get_int_argument(
		"cuba-vegas-increase-evals", default_integration_parameters.cuba.vegas.continuing_iterations
	);
	common_integration_parameters.cuba.verbose = args.get_int_argument("cuba-verbose", 0);

	integration_parameters lo_integration_parameters{common_integration_parameters};
	lo_integration_parameters.cuba.maximum_relative_error = args.get_double_argument(
		"cuba-max-rel-error-lo", common_integration_parameters.cuba.maximum_relative_error
	);

	integration_parameters nlo_integration_parameters{common_integration_parameters};
	nlo_integration_parameters.cuba.maximum_relative_error = args.get_double_argument(
		"cuba-max-rel-error-nlo", common_integration_parameters.cuba.maximum_relative_error
	);

	integration_parameters nnlo_integration_parameters{common_integration_parameters};
	nnlo_integration_parameters.cuba.maximum_relative_error = args.get_double_argument(
		"cuba-max-rel-error-nnlo", common_integration_parameters.cuba.maximum_relative_error
	);

	const std::unordered_map<perturbative_order, dm::integration_parameters> integration_parameters = {
		{perturbative_order::lo, lo_integration_parameters},
		{perturbative_order::nlo, nlo_integration_parameters},
		{perturbative_order::nnlo, nnlo_integration_parameters}
	};

	const bool show_only_channels = args.get_flag("show-only-channels");

	const unsigned int ff_member = args.get_unsigned_int_argument("ff-member", 0);

	// The kernel mode computes all interpolation kernels in one run, in place of the PDF and the FFs (see `vector_kernel_interface`),
	// with the decays of a single resonance. The kernels are those of piecewise-linear interpolation on the given knots, one for each knot.
	// They are integrated with plain quasi-Monte Carlo, with --cuba-max-evals points.
	const std::vector<double> pdf_knots = mode == run_mode::kernel ? args.get_required_double_list("pdf-knots") : std::vector<double>{};
	const std::vector<double> ff_knots = mode == run_mode::kernel ? args.get_required_double_list("ff-knots") : std::vector<double>{};
	const std::optional<particle> resonance_particle = mode == run_mode::kernel
		? std::make_optional(args.get_required_choice_argument("resonance-particle", maps::resonance_particles))
		: std::nullopt;

	// Read flavor options
	const flavor::flavor_vector flavors = args.get_flavor_list("flavors", {flavor::u, flavor::d, flavor::c, flavor::s, flavor::b});
	const flavor::flavor_vector charm_flavors = args.get_flavor_list("charm-flavors", {flavor::d, flavor::c, flavor::s, flavor::b});

	const flavor::flavor_vector initial_flavors = args.get_flavor_list("initial-flavors", std::ranges::to<flavor::flavor_vector>(flavor::all_flavors));
	const flavor::flavor_vector fragmentation_flavors = args.get_flavor_list("fragmentation-flavors", std::ranges::to<flavor::flavor_vector>(flavor::all_flavors));

	// Read PDF options. The PDF is used only in the direct mode: the kernel mode integrates the interpolation kernels in its place.
	const auto make_pdf = [&]() {
		const std::string pdf_set = args.get_required_argument("pdf");

		const unsigned int pdf_member = args.get_unsigned_int_argument("pdf-member", 0);

		lha_interface pdf(pdf_set, pdf_member);

		if (args.get_flag("isospin")) {
			const double Z = args.get_double_argument("Z", 1.0);
			const double A = args.get_double_argument("A", 1.0);

			pdf.enable_isospin(Z, A);
		}

		return pdf;
	};
	const std::optional<lha_interface<>> pdf = mode == run_mode::direct ? std::make_optional(make_pdf()) : std::nullopt;

	// Read process type
	const bool neutrino = args.get_flag("neutrino");
	const bool antineutrino = args.get_flag("antineutrino");

	if (neutrino == antineutrino) {
		throw std::runtime_error("Exactly one of --neutrino and --antineutrino is required");
	}
	const scattering_process process(
		neutrino ? scattering_process::process_type::neutrino : scattering_process::process_type::antineutrino, 
		constants::particles::proton, constants::particles::neutrino
	);

	// Read kinematical point
	const double x = args.get_required_double_argument("x");
	const double y = args.get_required_double_argument("y");
	const double E = args.get_required_double_argument("E");
	const trf_kinematics kinematics = trf_kinematics::y_E_beam(x, y, E, process.target.mass, process.projectile.mass);
	if (!kinematics.is_valid()) {
		throw std::runtime_error(std::format("The kinematical point --x {} --y {} --E {} is outside the physical region", x, y, E));
	}

	// Read quark masses. The charm mass defaults to that of the PDF, so in the kernel mode it must be given.
	const double strange_mass = args.get_double_argument("strange-mass", 0.0);
	const double down_mass = args.get_double_argument("down-mass", 0.0);
	const double up_mass = args.get_double_argument("up-mass", 0.0);
	const double top_mass = args.get_double_argument("top-mass", 0.0);
	const double bottom_mass = args.get_double_argument("bottom-mass", 0.0);

	const double charm_mass = pdf ? args.get_double_argument("charm-mass", pdf->quark_mass(flavor::c)) : args.get_required_double_argument("charm-mass");

	// Read scale multipliers
	const double renormalization_scale_multiplier = args.get_double_argument("renormalization", 1.0);
	const double factorization_scale_multiplier = args.get_double_argument("factorization", 1.0);
	const double fragmentation_scale_multiplier = args.get_double_argument("fragmentation", 1.0);

	// Read the scale dependence: mu^2 = Q^2 (Q2) or mu^2 = Q^2 + mc^2 (Q2+mc2, the default), see `scale` below
	const bool scale_adds_charm_mass = args.get_choice_argument("scale", maps::scale_dependences, true);
	const bool renormalization_bounded = !args.get_flag("renormalization-unbounded");
	const bool factorization_bounded = !args.get_flag("factorization-unbounded");
	const bool fragmentation_bounded = !args.get_flag("fragmentation-unbounded");

	// Read the minimum energy of the muons from the decays
	const double min_E = args.get_double_argument("minE", 5.0);

	// Construct scale dependence function: Q^2, plus mc^2 with --scale Q2+mc2, bounded from below by mc^2 if bound_from_below
	const auto scale = [scale_adds_charm_mass](const double charm_mass_, const bool bound_from_below) {
		const double mc2 = math::power<2>(charm_mass_);
		const double lower_bound = bound_from_below ? mc2 : 0.0;
		const double mass_term = scale_adds_charm_mass ? mc2 : 0.0;

		return scale_dependence::function([mass_term](const trf_kinematics &kinematics_) noexcept { return kinematics_.Q2 + mass_term; }, lower_bound);
	};

	// Read decay parametrization
	const std::optional<double> decay_N = args.get_double_argument("decay-N");
	const std::optional<double> decay_alpha = args.get_double_argument("decay-alpha");
	const std::optional<double> decay_beta = args.get_double_argument("decay-beta");
	const std::optional<double> decay_gamma = args.get_double_argument("decay-gamma");

	const decay_parametrization parametrization = [=]() {
		if (decay_N && decay_alpha && decay_beta && decay_gamma) {
			return decay_parametrization{*decay_N, *decay_alpha, *decay_beta, *decay_gamma};
		}
		if (decay_N || decay_alpha || decay_beta || decay_gamma) {
			throw std::runtime_error("The decay parametrization needs all of --decay-N, --decay-alpha, --decay-beta and --decay-gamma");
		}
		return decay_parametrization::default_fit();
	}();

	const auto grid_fragmentation = grids::create_grid_fragmentation({ff_member, ff_member, ff_member, 0});
	const std::array parton_masses{
		top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
	};

	const parton_configuration charm_parton_configuration{charm_flavors, parton_masses};
	const parton_configuration total_parton_configuration{flavors, parton_masses};

	const std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> dis_charm_channels = {
		{perturbative_order::lo, 	channels::read_dis_channels("lo", charm_parton_configuration, process, args)},
		{perturbative_order::nlo, 	channels::read_dis_channels("nlo", charm_parton_configuration, process, args)}
	};
	const std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> dis_channels = {
		{perturbative_order::lo, 	channels::read_dis_channels("lo", total_parton_configuration, process, args)},
		{perturbative_order::nlo, 	channels::read_dis_channels("nlo", total_parton_configuration, process, args)}
	};

	// The channels of the fragmenting flavor that is masked in the process are left out (see `grids::mask_fragmentation_flavors`),
	// unless --disable-fragmentation-channel-mask is given
	const flavor::flavor_vector sidis_fragmentation_flavors = args.get_flag("disable-fragmentation-channel-mask")
		? fragmentation_flavors
		: grids::mask_fragmentation_flavors(fragmentation_flavors, process);

	const std::unordered_map<perturbative_order, std::vector<sidis::partonic_channel>> sidis_channels = {
		{perturbative_order::lo, 	channels::filter_channels(channels::read_sidis_channels("lo", total_parton_configuration, process, args), initial_flavors, sidis_fragmentation_flavors)},
		{perturbative_order::nlo, 	channels::filter_channels(channels::read_sidis_channels("nlo", total_parton_configuration, process, args), initial_flavors, sidis_fragmentation_flavors)},
		{perturbative_order::nnlo, 	channels::filter_channels(channels::read_sidis_channels("nnlo", total_parton_configuration, process, args), initial_flavors, sidis_fragmentation_flavors)}
	};

	// SIDIS is computed up to the highest order that has channels. The orders below it are computed too, and those without channels stay zero.
	const perturbative_order order = highest_order(sidis_channels);

	// Read the thresholds of the heavy quarks, which determine the number of active flavors (NF) in the NNLO coefficients.
	// These are distinct from the quark masses above, which set the final-state masses of the channels.
	// By default they are the quark masses of the PDF. The kernel mode has no PDF, so there they must be given if there are NNLO channels
	// (NF appears only in the NNLO coefficients), and they should be the quark masses of the PDF sets that the grids are convolved with.
	const auto read_threshold = [&](const std::string name, const flavor::type flavor) {
		if (pdf) {
			return args.get_double_argument(name, pdf->quark_mass(flavor));
		}
		return order >= perturbative_order::nnlo ? args.get_required_double_argument(name) : args.get_double_argument(name, 0.0);
	};
	const heavy_quark_thresholds thresholds{
		read_threshold("charm-threshold", flavor::c), read_threshold("bottom-threshold", flavor::b), read_threshold("top-threshold", flavor::t)
	};

	// Read the coefficient implementation of each perturbative order. The NNLO coefficient functions have only the gc implementation (checked by `sidis::instance`).
	const std::unordered_map<perturbative_order, coefficient_implementation_t> coefficient_implementation = {
		{perturbative_order::lo, 	args.get_choice_argument("coefficient-implementation-lo", maps::coefficient_implementations, coefficient_implementation_t::hw)},
		{perturbative_order::nlo, 	args.get_choice_argument("coefficient-implementation-nlo", maps::coefficient_implementations, coefficient_implementation_t::hw)},
		{perturbative_order::nnlo, 	args.get_choice_argument("coefficient-implementation-nnlo", maps::coefficient_implementations, coefficient_implementation_t::gc)}
	};

	const cross_section_observable observable;

	// All arguments have been read
	args.check_arguments();

	if (observable_name == "dis.differential.charm" || observable_name == "dis.differential.total") {
		if (mode != run_mode::direct) {
			throw std::runtime_error("Only direct mode is supported with DIS");
		}

		const bool charm_dis = observable_name == "dis.differential.charm";

		const std::unordered_map<perturbative_order, std::vector<dis::partonic_channel>> channels = charm_dis ? dis_charm_channels : dis_channels;

		// The DIS coefficient functions have only the hw implementation
		const std::unordered_map<perturbative_order, coefficient_implementation_t> dis_coefficient_implementation = {
			{perturbative_order::lo, coefficient_implementation_t::hw},
			{perturbative_order::nlo, coefficient_implementation_t::hw}
		};

		check_channels(channels);
		print_channels(channels, dis_coefficient_implementation);
		if (show_only_channels) { return EXIT_SUCCESS; }

		const auto renormalization_scale_function = renormalization_scale_multiplier * scale(charm_mass, renormalization_bounded);
		const auto factorization_scale_function = factorization_scale_multiplier * scale(charm_mass, factorization_bounded);

		// DIS is computed up to NLO, the highest order of its coefficient functions
		const dis::instance dis{
			ckm_matrix::global_fit(),
			*pdf,
			integration_parameters,
			process,
			renormalization_scale_function, factorization_scale_function,
			true, observable,
			perturbative_order::nlo,
			0.0,
			channels
		};

		const dis::result result = dis.differential_cross_section_xy(kinematics);

		const double renormalization_scale = renormalization_scale_function(kinematics);
		const double factorization_scale = factorization_scale_function(kinematics);

		const double alpha_s = dis.compute_alpha_s(kinematics);

		// The blocks of the SIDIS output (the totals, a row per physical channel, and the sums of the rows), without the fragmentation scale and flavor and up to NLO
		write_output(output, output_file, [&](std::ostream &stream) {
			stream << "===" << std::endl;

			stream << pdf->get_set_name() << ", ";
			stream << pdf->get_set_member() << ", ";
			
			stream << kinematics.x << ", ";
			stream << kinematics.y << ", ";
			stream << kinematics.E_beam << ", ";
			
			stream << kinematics.Q2 << ", ";
			stream << renormalization_scale << ", ";
			stream << factorization_scale << ", ";
			stream << alpha_s << ", ";
			
			stream << result.value.lo() << ", ";
			stream << result.value.nlo() << ", ";
			stream << result.value.nlo_correction_with_alpha_s() << std::endl;

			double lo_sum = 0.0;
			double nlo_sum = 0.0;

			stream << "===" << std::endl;
			for (const auto &[physical_channel, value] : std::views::zip(result.physical_channels, result.channel_decomposition)) {
				stream << std::to_underlying(physical_channel.initial_flavor) << ", ";
				stream << value.lo() << ", ";
				stream << value.nlo_correction() << std::endl;
				lo_sum += value.lo_correction();
				nlo_sum += value.nlo_correction();
			}
			stream << "===" << std::endl;

			stream << lo_sum << std::endl;
			stream << nlo_sum << std::endl;

			stream << "===" << std::endl;
		});
	} else if (observable_name == "sidis.differential") {
		check_channels(sidis_channels);
		print_channels(sidis_channels, coefficient_implementation);
		if (show_only_channels) { return EXIT_SUCCESS; }

		const auto renormalization_scale_function = renormalization_scale_multiplier * scale(charm_mass, renormalization_bounded);
		const auto factorization_scale_function = factorization_scale_multiplier * scale(charm_mass, factorization_bounded);
		const auto fragmentation_scale_function = fragmentation_scale_multiplier * scale(charm_mass, fragmentation_bounded);

		// The instance with the given PDF and FFs, the FFs in the same order as the decays
		const auto make_instance = [&](const auto &pdf_interface, const auto &ffs, const auto &decays) {
			return sidis::instance{
				ckm_matrix::global_fit(),
				pdf_interface, ffs, decays,
				integration_parameters,
				process,
				renormalization_scale_function, factorization_scale_function, fragmentation_scale_function,
				true, observable,
				order,
				0.0,
				sidis_channels,
				coefficient_implementation,
				thresholds
			};
		};

		switch (mode) {
			case run_mode::direct: {
				const auto [ffs, decays] = grid_fragmentation(min_E, constants::particles::massless_muon, parametrization, process);
				const auto sidis_instance = make_instance(*pdf, ffs, decays);

				const sidis::result result = sidis_instance.lepton_pair_cross_section_xy(kinematics);

				const double renormalization_scale = renormalization_scale_function(kinematics);
				const double factorization_scale = factorization_scale_function(kinematics);
				const double fragmentation_scale = fragmentation_scale_function(kinematics);
		
				const double alpha_s = sidis_instance.compute_alpha_s(kinematics);
		
				write_output(output, output_file, [&](std::ostream &stream) {
					stream << "===" << std::endl;

					stream << pdf->get_set_name() << ", ";
					stream << pdf->get_set_member() << ", ";
					
					stream << kinematics.x << ", ";
					stream << kinematics.y << ", ";
					stream << kinematics.E_beam << ", ";
					
					stream << kinematics.Q2 << ", ";
					stream << renormalization_scale << ", ";
					stream << factorization_scale << ", ";
					stream << fragmentation_scale << ", ";
					stream << alpha_s << ", ";
					
					stream << result.value.lo() << ", ";
					stream << result.value.nlo() << ", ";
					stream << result.value.nnlo() << ", ";
					stream << result.value.nlo_correction_with_alpha_s() << ", ";
					stream << result.value.nnlo_correction_with_alpha_s() << std::endl;

					double lo_sum = 0.0;
					double nlo_sum = 0.0;
					double nnlo_sum = 0.0;

					stream << "===" << std::endl;
					for (const auto &[physical_channel_index, physical_channel] : std::views::enumerate(result.physical_channels)) {
						// The PDF and FF values have a single row each (see `sidis::parton_values`)
						const perturbative_quantity &value = result.channel_decomposition[0, 0, physical_channel_index];
						stream << std::to_underlying(physical_channel.initial_flavor) << ", ";
						stream << std::to_underlying(physical_channel.fragmentation_flavor) << ", ";
						stream << value.lo() << ", ";
						stream << value.nlo_correction() << ", ";
						stream << value.nnlo_correction() << std::endl;
						lo_sum += value.lo_correction();
						nlo_sum += value.nlo_correction();
						nnlo_sum += value.nnlo_correction();
					}
					stream << "===" << std::endl;

					stream << lo_sum << std::endl;
					stream << nlo_sum << std::endl;
					stream << nnlo_sum << std::endl;

					stream << "===" << std::endl;
				});

				break;
			}

			case run_mode::kernel: {
				// The interpolation kernels take the place of the PDF and the FFs (see `vector_kernel_interface`), with the decays of a single resonance.
				// The kernels don't depend on the FF sets, so the per-flavor multipliers of the FF sets (see `grids::create_grid_fragmentation`) are not applied.
				const auto resonance_decays = grids::create_resonance_decays(*resonance_particle);
				const auto decays = resonance_decays(min_E, constants::particles::massless_muon, parametrization, process);
				const vector_kernel_interface ff_kernels{linear_interpolation_basis{ff_knots}};
				const auto sidis_instance = make_instance(vector_kernel_interface{linear_interpolation_basis{pdf_knots}}, std::vector(decays.size(), ff_kernels), decays);

				const sidis::result grid = sidis_instance.lepton_pair_cross_section_xy(kinematics);

				write_output(output, output_file, [&](std::ostream &stream) {
					// The kinematics and the scales
					stream << "===" << std::endl;
					stream << kinematics.x << ", ";
					stream << kinematics.y << ", ";
					stream << kinematics.E_beam << ", ";
					stream << kinematics.Q2 << ", ";
					stream << renormalization_scale_function(kinematics) << ", ";
					stream << factorization_scale_function(kinematics) << ", ";
					stream << fragmentation_scale_function(kinematics) << std::endl;

					// One row per pair of kernels and physical channel: the kernels, the flavors, and the values of the orders, as in the channel decomposition
					stream << "===" << std::endl;
					for (kernel_index pdf_kernel = 0; pdf_kernel < grid.channel_decomposition.extent(0); ++pdf_kernel) {
						for (kernel_index ff_kernel = 0; ff_kernel < grid.channel_decomposition.extent(1); ++ff_kernel) {
							for (const auto &[physical_channel_index, physical_channel] : std::views::enumerate(grid.physical_channels)) {
								const perturbative_quantity &value = grid.channel_decomposition[pdf_kernel, ff_kernel, physical_channel_index];
								stream << pdf_kernel << ", ";
								stream << ff_kernel << ", ";
								stream << std::to_underlying(physical_channel.initial_flavor) << ", ";
								stream << std::to_underlying(physical_channel.fragmentation_flavor) << ", ";
								stream << value.lo() << ", ";
								stream << value.nlo_correction() << ", ";
								stream << value.nnlo_correction() << std::endl;
							}
						}
					}
					stream << "===" << std::endl;
				});

				break;
			}

			std::unreachable();
		}
	} else {
		throw std::runtime_error(std::format("Unknown observable '{}'", observable_name));
	}

	return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
	// Errors, such as missing or invalid arguments, are printed as messages rather than terminating with an uncaught exception
	try {
		return run(arg_parser(argc, argv));
	} catch (const std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}
