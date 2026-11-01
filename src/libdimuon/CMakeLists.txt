# Create the header-only library
add_library(libdimuon INTERFACE)
# Alias `libdimuon` to the conventional `dimuon::dimuon`
add_library(dimuon::dimuon ALIAS libdimuon)
set_target_properties(libdimuon PROPERTIES EXPORT_NAME dimuon)
# Require c++23
target_compile_features(libdimuon INTERFACE cxx_std_23)

# Define file sets for headers and fragment files (which do not compile on their own - see `dimuon/sidis/coefficients/gc`)
target_sources(libdimuon INTERFACE FILE_SET HEADERS FILE_SET fragments TYPE HEADERS)
# Include only the headers (and not the fragment files) into the header verification target
set_target_properties(libdimuon PROPERTIES INTERFACE_HEADER_SETS_TO_VERIFY HEADERS)
# Include the library subdirectories, which actually define the headers and fragment files
add_subdirectory(dimuon/common)
add_subdirectory(dimuon/coupling)
add_subdirectory(dimuon/decay)
add_subdirectory(dimuon/dis)
add_subdirectory(dimuon/integration)
add_subdirectory(dimuon/interpolation)
add_subdirectory(dimuon/pdf)
add_subdirectory(dimuon/sidis)
add_subdirectory(dimuon/utility)

# Include third-party code as system headers so that 
# the warning flags used in this project are not applied to them
target_include_directories(
	libdimuon SYSTEM INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/dimuon/external>
)

# Add dependencies
find_package(GSL REQUIRED)
find_package(LHAPDF REQUIRED)
find_package(Cuba REQUIRED)
target_link_libraries(libdimuon INTERFACE GSL::gsl LHAPDF::LHAPDF Cuba::Cuba)

# Set the optimization flags
# This is a header-only library: it's not compiled on its own, but only when linked with another target
# that is actually compiled (like the dimuon executable). Thus any optimization flags have to come from
# that other target. This defines a variable that contains the recommended flags that will be exported
# for easy access by users.
set(optimization_flags -fno-math-errno -fno-trapping-math -fassociative-math -fno-signed-zeros)
# Set `-march=native` unless `DIMUON_NATIVE` has been set to `OFF`. This flag, among other things, selects
# the appropriate vector instruction set and can therefore have an impact on performance. However, it also
# makes the build significantly less portable.
if (DIMUON_NATIVE)
	list(APPEND optimization_flags -march=native)
endif()
# Create an empty target for storing the optimization flags
add_library(dimuon_optimization_flags INTERFACE)
# Alias the empty target to `dimuon::optimization_flags`
add_library(dimuon::optimization_flags ALIAS dimuon_optimization_flags)
set_target_properties(dimuon_optimization_flags PROPERTIES EXPORT_NAME optimization_flags)
# Attach the flags to the empty target via `INTERFACE`, so that any other target 
# linking against `dimuon::optimization_flags` will inherit those same options.
# The flags are included conditionally, so that they are not applied in a debug build
target_compile_options(dimuon_optimization_flags INTERFACE "$<$<NOT:$<CONFIG:Debug>>:${optimization_flags}>")

# Install (i.e. copy) the headers.
# Invoked with `cmake --install <build directory> --prefix <directory>`. 
# Headers are installed to <directory>/include/dimuon. 
# Third-party code and licenses are installed to <directory>/include/dimuon/external.
install(
	TARGETS libdimuon dimuon_optimization_flags EXPORT dimuonTargets
	FILE_SET HEADERS FILE_SET fragments
)
install(DIRECTORY external DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dimuon)

# Generate the `dimuonTargets.cmake` package that allows users to write `find_package(dimuon)`
include(CMakePackageConfigHelpers)
write_basic_package_version_file(dimuonConfigVersion.cmake COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(EXPORT dimuonTargets NAMESPACE dimuon:: DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/dimuon)
install(
	FILES
		${PROJECT_SOURCE_DIR}/cmake/dimuonConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/dimuonConfigVersion.cmake
		${PROJECT_SOURCE_DIR}/cmake/FindLHAPDF.cmake
		${PROJECT_SOURCE_DIR}/cmake/FindCuba.cmake
	DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/dimuon
)
