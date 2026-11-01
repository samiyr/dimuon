add_subdirectory(gc)
add_subdirectory(hw)

target_sources(libdimuon
	INTERFACE FILE_SET HEADERS FILES
		coefficient_function.hpp
		evaluation_parameters.hpp
		f1_hw.hpp
		f2_hw.hpp
		f3_hw.hpp
		f1_gc.hpp
		f2_gc.hpp
		f3_gc.hpp
		evaluation_cache_gc.hpp
)
