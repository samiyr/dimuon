#[=======================================================================[.rst:
FindLHAPDF
----------

Finds the LHAPDF library (https://www.lhapdf.org).

The library is searched for below ``LHAPDF_ROOT`` or ``LHAPDF_ROOT_DIR`` (CMake or environment variables), then below
the prefix that ``lhapdf-config`` reports if that program is on the ``PATH``, and then in the default locations.

Imported target: ``LHAPDF::LHAPDF``

Result variables: ``LHAPDF_FOUND``, ``LHAPDF_VERSION``

Cache variables: ``LHAPDF_INCLUDE_DIR``, ``LHAPDF_LIBRARY``
#]=======================================================================]

find_program(LHAPDF_CONFIG_EXECUTABLE lhapdf-config)
if (LHAPDF_CONFIG_EXECUTABLE)
	execute_process(
		COMMAND ${LHAPDF_CONFIG_EXECUTABLE} --prefix
		OUTPUT_VARIABLE LHAPDF_config_prefix OUTPUT_STRIP_TRAILING_WHITESPACE
	)
endif()

set(LHAPDF_hints ${LHAPDF_ROOT_DIR} ENV LHAPDF_ROOT_DIR ${LHAPDF_config_prefix})
find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h HINTS ${LHAPDF_hints} PATH_SUFFIXES include)
find_library(LHAPDF_LIBRARY LHAPDF HINTS ${LHAPDF_hints} PATH_SUFFIXES lib lib64)
mark_as_advanced(LHAPDF_CONFIG_EXECUTABLE LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY)

if (LHAPDF_INCLUDE_DIR)
	file(STRINGS ${LHAPDF_INCLUDE_DIR}/LHAPDF/Version.h LHAPDF_version_line REGEX "#define LHAPDF_VERSION \"")
	string(REGEX MATCH "[0-9.]+" LHAPDF_VERSION "${LHAPDF_version_line}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	LHAPDF
	REQUIRED_VARS LHAPDF_LIBRARY LHAPDF_INCLUDE_DIR
	VERSION_VAR LHAPDF_VERSION
)

if (LHAPDF_FOUND AND NOT TARGET LHAPDF::LHAPDF)
	add_library(LHAPDF::LHAPDF UNKNOWN IMPORTED)
	set_target_properties(
		LHAPDF::LHAPDF PROPERTIES
		IMPORTED_LOCATION ${LHAPDF_LIBRARY}
		INTERFACE_INCLUDE_DIRECTORIES ${LHAPDF_INCLUDE_DIR}
	)
endif()
