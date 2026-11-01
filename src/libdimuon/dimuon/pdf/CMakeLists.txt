add_subdirectory(extrapolator)
add_subdirectory(interface)

target_sources(libdimuon
	INTERFACE FILE_SET HEADERS FILES
		is_pdf_interface.hpp
)
