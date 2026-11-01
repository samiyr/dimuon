from enum import Enum

class Observable(Enum):
	"""
	Enum describing the observable type, used by the CC program
	"""

	charm_dis = "dis.differential.charm"
	sidis = "sidis.differential"

	def __str__(self):
		return self.value

class RunMode(Enum):
	"""
	Enum describing how the CC program computes the observable (its --run-mode argument)
	"""

	direct = "direct"
	"""
	"Normal" run mode, where the PDF and FFs are supplied, and the program computes full cross section values
	"""

	kernel = "kernel"
	"""
	Computes interpolation tables
	"""

	def __str__(self):
		return self.value
	
class Mode(Enum):
	"""
	Enum describing the variations that dimuon_cross_sections.py computes and how their outputs are split (the --mode argument)
	"""

	normal = "normal"
	"""
	In the normal mode, PDF error sets are split into separate files
	"""

	masses = "masses"
	"""
	In the masses mode, different charm masses are split into separate folders
	"""

	scales7p = "scales7p"
	"""
	In the scales mode, scale variations are enabled and they are split into separate files.
	In this mode, the fragmentation scale is kept at 1, resulting in 7 scale variations. Can be used for DIS and SIDIS.
	If multiple PDF set members are included, they will be ignored and only the last PDF member data will be kept.
	"""

	scales17p = "scales17p"
	"""
	In the scales mode, scale variations are enabled and they are split into separate files.
	In this mode, the fragmentation scale is also varied, resulting in 17 scale variations. Can be used for SIDIS only.
	If multiple PDF set members are included, they will be ignored and only the last PDF member data will be kept.
	"""

	decay_function_errors = "decay-function-errors"
	"""
	In the decay-function error mode, a set of different decay-function parametrizations is used and split into separate files.
	This enables the propagation of decay-function errors using replicas.
	"""

	def __str__(self):
		return self.value
	
class ScaleDependence(Enum):
	"""
	Enum describing scale dependence options
	"""

	Q2 = "Q2"
	"""
	mu^2 = max(m_c^2, k Q^2), where k is the scale multiplier (usually 0.25, 1.0, or 4.0)
	"""

	Q2mc2 = "Q2+mc2"
	"""
	mu^2 = max(m_c^2, k(Q^2 + m_c^2)), where k is the scale multiplier (usually 0.25, 1.0, or 4.0)
	"""

	def __str__(self):
		return self.value

class Channels(Enum):
	dis_zm = "dis-zm"
	dis_zm_q = "dis-zm-q"
	dis_zm_g = "dis-zm-g"

	dis_sr = "dis-sr"
	dis_sr_q = "dis-sr-q"
	dis_sr_g = "dis-sr-g"
	
	dis_gm = "dis-gm"
	dis_gm_neutrino = "dis-gm-neutrino"
	dis_gm_antineutrino = "dis-gm-antineutrino"

	dis_gm_unmodified = "dis-gm-unmod"
	dis_gm_neutrino_unmodified = "dis-gm-neutrino-unmod"
	dis_gm_antineutrino_unmodified = "dis-gm-antineutrino-unmod"

	dis_gm_fixed_order = "dis-gm-fo"
	dis_gm_neutrino_fixed_order = "dis-gm-neutrino-fo"
	dis_gm_antineutrino_fixed_order = "dis-gm-antineutrino-fo"

	sidis_zm = "sidis-zm"
	sidis_sr = "sidis-sr"

	sidis_gm = "sidis-gm"

	def __str__(self):
		return self.value

class CoefficientImplementation(Enum):
	hw = "hw"
	gc = "gc"

	def __str__(self):
		return self.value

class IntegrationMethod(Enum):
	"""
	Enum describing the integration method of the CC program in the direct mode (its --integration-method argument).
	The kernel mode integrates with plain quasi-Monte Carlo.
	"""

	suave = "suave"
	"""
	Cuba's Suave
	"""

	vegas = "vegas"
	"""
	Cuba's Vegas
	"""

	qmc = "qmc"
	"""
	Plain quasi-Monte Carlo with Sobol points, as many as --cuba-max-evals
	"""

	def __str__(self):
		return self.value

class ResonanceParticle(Enum):
	D0 = "D0"
	Dp = "Dp"
	Ds = "Ds"
	LambdaC = "LambdaC"

	def __str__(self):
		return self.value