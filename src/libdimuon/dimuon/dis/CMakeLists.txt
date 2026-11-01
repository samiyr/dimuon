add_subdirectory(coefficients)

target_sources(libdimuon
	INTERFACE FILE_SET HEADERS FILES
		channel_factory.hpp
		coefficient_channel.hpp
		evaluator.hpp
		instance.hpp
		parameters.hpp
		parton_values.hpp
		partonic_channel.hpp
		physical_channel.hpp
		result.hpp
		vectorization.hpp
)
