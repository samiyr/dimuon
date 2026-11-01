import numpy as np
from .arg_modes import Observable, Mode, Channels

nutev_neutrino_x = [0.021, 0.057, 0.101, 0.167, 0.326]
nutev_neutrino_y = [0.324, 0.558, 0.771]
nutev_neutrino_E = [88.29, 174.29, 247.0]

nutev_antineutrino_x = [0.015, 0.042, 0.072, 0.114, 0.207]
nutev_antineutrino_y = [0.349, 0.579, 0.776]
nutev_antineutrino_E = [77.9, 143.7, 226.8]

ccfr_neutrino_x = [0.023, 0.057, 0.100, 0.167, 0.336]
ccfr_neutrino_y = [0.32, 0.57, 0.795]
ccfr_neutrino_E = [109.46, 209.89, 332.7]

ccfr_antineutrino_x = [0.018, 0.041, 0.069, 0.111, 0.210]
ccfr_antineutrino_y = [0.355, 0.596, 0.802]
ccfr_antineutrino_E = [87.40, 160.52, 265.76]

x_bins_1 = np.linspace(0.01, 0.1, 91)
x_bins_2 = np.linspace(0.11, 0.4, 30)
x_bins = np.concatenate([x_bins_1, x_bins_2])

def get_df_columns(observable:Observable):
	"""
	Returns a list of dataframe column names used before splicing the single output dataframe into several smaller ones.

	Args:
		observable: observable that the CC program computes
	"""

	dis_differential_columns = [
		"pdf_set", "pdf_member", 
		"x", "y", "E",
		"Q2", "renormalization_scale", "factorization_scale", "alpha_s",
		"LO", "NLO", "NLO_correction",
		"category", "charm_mass", "scales", "parametrization"
	]

	sidis_differential_columns = [
		"pdf_set", "pdf_member", 
		"x", "y", "E",
		"Q2", "renormalization_scale", "factorization_scale", "fragmentation_scale", "alpha_s",
		"LO", "NLO", "NNLO", "NLO_correction", "NNLO_correction",
		"category", "charm_mass", "scales", "parametrization"
	]

	columns = {
		Observable.charm_dis: dis_differential_columns,
		Observable.sidis: sidis_differential_columns
	}

	return columns[observable]

def get_output_columns(args):
	"""
	Returns a list of dataframe column names that appear in the exported CSV files.

	Args:
		args: command line arguments from argparse
	"""
	if args.observable == Observable.charm_dis:
		return [
			"x", "y", "E", 
			"Q2", "renormalization_scale", "factorization_scale", "alpha_s",
			"LO", "NLO", "NLO_correction",
		]
	elif args.observable == Observable.sidis:
		return [
			"x", "y", "E",
			"Q2", "renormalization_scale", "factorization_scale", "fragmentation_scale", "alpha_s",
			"LO", "NLO", "NNLO", "NLO_correction", "NNLO_correction"
		]

def get_data_folder(args):
	"""
	Constructs the data folder path based on the given command line arguments.

	Args:
		args: command line arguments from argparse
	"""
	starts = {
		Observable.charm_dis: "direct-mode/charm_production",
		Observable.sidis: "direct-mode/dimuon_production"
	}
	start = starts[args.observable]

	observable_folders = {
		Mode.normal: "error_sets",
		Mode.masses: "charm_masses",
		Mode.scales7p: "scales7p",
		Mode.scales17p: "scales17p",
		Mode.decay_function_errors: "decay_function_error_sets"
	}
	observable_folder = observable_folders[args.mode]

	return f"{start}/{observable_folder}"

channel_presets = {
	Channels.dis_zm: [
		"--lo-all-m"	, "zm",
		"--lo-all-M"	, "zm",
		"--nlo-all-m"	, "zm",
		"--nlo-all-M"	, "zm"
	],
	Channels.dis_zm_q: [
		"--lo-q-m"		, "zm",
		"--lo-qbar-m"	, "zm",
		"--nlo-q-m"		, "zm",
		"--nlo-qbar-m"	, "zm"
	],
	Channels.dis_zm_g: [
		"--nlo-g-m"		, "zm",
		"--nlo-gbar-m"	, "zm"
	],
	Channels.dis_sr: [
		"--lo-all-m"	, "sr",
		"--lo-all-M"	, "sr",
		"--nlo-all-m"	, "sr",
		"--nlo-all-M"	, "sr"
	],
	Channels.dis_sr_q: [
		"--lo-q-m"		, "sr",
		"--lo-q-M"		, "sr",
		"--lo-qbar-m"	, "sr",
		"--lo-qbar-M"	, "sr",
		"--nlo-q-m"		, "sr",
		"--nlo-q-M"		, "sr",
		"--nlo-qbar-m"	, "sr",
		"--nlo-qbar-M"	, "sr"
	],
	Channels.dis_sr_g: [
		"--nlo-g-m"		, "sr",
		"--nlo-g-M"		, "sr",
		"--nlo-gbar-m"	, "sr",
		"--nlo-gbar-M"	, "sr"
	],
	Channels.dis_gm_neutrino: [
		"--lo-q-m"		, "sr",
		"--lo-q-M"		, "gm",
		"--lo-qbar-m"	, "sr*",
		"--lo-qbar-M"	, "gm*",
		"--nlo-q-m"		, "sr",
		"--nlo-q-M"		, "gm",
		"--nlo-qbar-m"	, "sr*",
		"--nlo-qbar-M"	, "gm*",
		"--nlo-g-m"		, "sr*",
		"--nlo-g-M"		, "gm*",
		"--nlo-gbar-m"	, "sr*",
		"--nlo-gbar-M"	, "gm*",
	],
	Channels.dis_gm_antineutrino: [
		"--lo-q-m"		, "sr*",
		"--lo-q-M"		, "gm*",
		"--lo-qbar-m"	, "sr",
		"--lo-qbar-M"	, "gm",
		"--nlo-q-m"		, "sr*",
		"--nlo-q-M"		, "gm*",
		"--nlo-qbar-m"	, "sr",
		"--nlo-qbar-M"	, "gm",
		"--nlo-g-m"		, "sr*",
		"--nlo-g-M"		, "gm*",
		"--nlo-gbar-m"	, "sr*",
		"--nlo-gbar-M"	, "gm*",
	],
	Channels.dis_gm_neutrino_unmodified: [
		"--lo-all-m"	, "sr",
		"--lo-all-M"	, "gm",
		"--nlo-all-m"	, "sr",
		"--nlo-all-M"	, "gm"
	],
	Channels.dis_gm_antineutrino_unmodified: [
		"--lo-all-m"	, "sr",
		"--lo-all-M"	, "gm",
		"--nlo-all-m"	, "sr",
		"--nlo-all-M"	, "gm"
	],
	Channels.dis_gm_neutrino_fixed_order: [
		"--lo-q-m"		, "sr",
		"--lo-q-M"		, "gm",
		"--nlo-q-m"		, "sr",
		"--nlo-q-M"		, "gm",
		"--nlo-g-m"		, "sr",
		"--nlo-g-M"		, "gm",
		"--nlo-gbar-m"	, "sr",
		"--nlo-gbar-M"	, "gm",
	],
	Channels.dis_gm_antineutrino_fixed_order: [
		"--lo-qbar-m"	, "sr",
		"--lo-qbar-M"	, "gm",
		"--nlo-qbar-m"	, "sr",
		"--nlo-qbar-M"	, "gm",
		"--nlo-g-m"		, "sr",
		"--nlo-g-M"		, "gm",
		"--nlo-gbar-m"	, "sr",
		"--nlo-gbar-M"	, "gm",
	],

	Channels.sidis_zm: [
		"--lo-all-m"	, "zm",
		"--lo-all-M"	, "zm",
		"--nlo-all-m"	, "zm",
		"--nlo-all-M"	, "zm"
	],
	Channels.sidis_sr: [
		"--lo-all-m"	, "sr",
		"--lo-all-M"	, "sr",
		"--nlo-all-m"	, "sr",
		"--nlo-all-M"	, "sr",
		"--nnlo-all-m"	, "sr",
		"--nnlo-all-M"	, "sr"
	],
	Channels.sidis_gm: [
		"--lo-all-m"	, "sr",
		"--lo-all-M"	, "gm",
		"--nlo-all-m"	, "sr",
		"--nlo-all-M"	, "gm",
		"--nnlo-all-m"	, "sr",
		"--nnlo-all-M"	, "sr"
	]
}

def get_channels(neutrino:bool, args):
	if args.channels:
		channel_preset = args.channels

		if channel_preset == Channels.dis_gm:
			return channel_presets[Channels.dis_gm_neutrino] if neutrino else channel_presets[Channels.dis_gm_antineutrino]
		elif channel_preset == Channels.dis_gm_unmodified:
			return channel_presets[Channels.dis_gm_neutrino_unmodified] if neutrino else channel_presets[Channels.dis_gm_antineutrino_unmodified]
		elif channel_preset == Channels.dis_gm_fixed_order:
			return channel_presets[Channels.dis_gm_neutrino_fixed_order] if neutrino else channel_presets[Channels.dis_gm_antineutrino_fixed_order]
		else:
			return channel_presets[channel_preset]
	else:
		return args.channel_list