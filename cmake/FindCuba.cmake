#[=======================================================================[.rst:
FindCuba
--------

Finds the Cuba library for multidimensional numerical integration (https://feynarts.de/cuba).

The library is searched for below ``Cuba_ROOT`` or ``CUBA_ROOT`` (CMake or environment variables), and then in the
default locations.

Imported target: ``Cuba::Cuba``

Result variables: ``Cuba_FOUND``

Cache variables: ``Cuba_INCLUDE_DIR``, ``Cuba_LIBRARY``
#]=======================================================================]

find_path(Cuba_INCLUDE_DIR cuba.h)
find_library(Cuba_LIBRARY cuba)
mark_as_advanced(Cuba_INCLUDE_DIR Cuba_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cuba REQUIRED_VARS Cuba_LIBRARY Cuba_INCLUDE_DIR)

if (Cuba_FOUND AND NOT TARGET Cuba::Cuba)
	add_library(Cuba::Cuba UNKNOWN IMPORTED)
	set_target_properties(
		Cuba::Cuba PROPERTIES
		IMPORTED_LOCATION ${Cuba_LIBRARY}
		INTERFACE_INCLUDE_DIRECTORIES ${Cuba_INCLUDE_DIR}
	)
endif()
