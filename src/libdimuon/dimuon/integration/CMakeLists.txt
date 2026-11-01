target_sources(libdimuon
	INTERFACE FILE_SET HEADERS FILES
		integration_parameters.hpp
		integration_result.hpp
		integrator_cquad.hpp
		integrator_mc.hpp
		integrator_nd.hpp
		is_1d_integrand_function.hpp
		is_1d_integrator.hpp
)
