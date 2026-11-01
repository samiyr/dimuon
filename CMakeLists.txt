cmake_minimum_required(VERSION 3.31)

# Define the project
project(
	dimuon
	VERSION 1.0.0
	DESCRIPTION "Dimuon production in neutrino deep inelastic scattering"
	LANGUAGES CXX
)

# Support for standard GNU installation directories
include(GNUInstallDirs)
# Supporting files for finding dependencies that aren't otherwise findable by CMake
list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Enable testing here so that ctest can find them.
# The tests themselves are added in `src` when configured with `-DBUILD_TESTS=ON`.
enable_testing()

# Add the library, the driver, the tools, and the tests
add_subdirectory(src)

# Install the decay grids to <directory>/share/dimuon/decaygrids
install(DIRECTORY decaygrids DESTINATION ${CMAKE_INSTALL_DATADIR}/dimuon)
# Install the examples to <directory>/share/dimuon/examples
install(DIRECTORY examples DESTINATION ${CMAKE_INSTALL_DATADIR}/dimuon)
# Install the runcards to <directory>/share/dimuon/runcards
install(DIRECTORY runcards DESTINATION ${CMAKE_INSTALL_DATADIR}/dimuon)
# Install the documentation to <directory>/share/doc/dimuon
install(FILES numerical_evaluation.md DESTINATION ${CMAKE_INSTALL_DOCDIR})