import datetime

from . import constants

def write(file, args):
	"""
	Writes a timestamp, all arguments of the run, and the channels passed to the CC program to an output file as comment lines.

	Args:
		file: Output file opened in text mode
		args: Arguments parsed by arg_parser
	"""
	file.write(f"# timestamp: {datetime.datetime.now().isoformat()}\n")

	for key, value in vars(args).items():
		value_str = str(value)
		if isinstance(value, list):
			value_str = ';'.join(map(str, value))
		file.write(f"# {key}: {value_str}\n")

	if args.neutrino:
		file.write("# neutrino channels:\n")
		channels = constants.get_channels(True, args)
		for channel, spec in zip(*[iter(channels)] * 2):
			file.write(f"# {channel}: {spec}\n")
	if args.antineutrino:
		file.write("# antineutrino channels:\n")
		channels = constants.get_channels(False, args)
		for channel, spec in zip(*[iter(channels)] * 2):
			file.write(f"# {channel}: {spec}\n")
