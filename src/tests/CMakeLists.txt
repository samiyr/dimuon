# The tests of the library (see README.md)

# GoogleTest: an installed copy if CMake finds one, otherwise downloaded
include(FetchContent)
FetchContent_Declare(
	googletest
	URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
	URL_HASH SHA256=65fab701d9829d38cb77c14acdc431d2108bfdbf8979e40eb8ae567edf10b27c
	FIND_PACKAGE_ARGS NAMES GTest
)
# GoogleTest is not installed with dimuon, and is compiled without the warning flags of this project
set(INSTALL_GTEST OFF)
get_directory_property(project_compile_options COMPILE_OPTIONS)
set_directory_properties(PROPERTIES COMPILE_OPTIONS "")
FetchContent_MakeAvailable(googletest)
set_directory_properties(PROPERTIES COMPILE_OPTIONS "${project_compile_options}")
include(GoogleTest)

# Fast tests of individual components
add_executable(unit_tests
	unit/decay_tests.cpp
	unit/dis_channel_tests.cpp
	unit/flavor_tests.cpp
	unit/instance_tests.cpp
	unit/integration_tests.cpp
	unit/interpolation_tests.cpp
	unit/kinematics_tests.cpp
	unit/pdf_tests.cpp
	unit/sidis_channel_tests.cpp
	unit/utility_tests.cpp
)
# Tests of physical properties of the cross sections: symmetries, limits, and the relation between SIDIS and DIS
add_executable(physics_tests
	physics/dis_conjugation_tests.cpp
	physics/dis_mass_scheme_tests.cpp
	physics/sidis_conjugation_tests.cpp
	physics/sidis_kernel_grid_tests.cpp
	physics/sidis_mass_scheme_tests.cpp
	physics/sum_rule_tests.cpp
)
# Comparisons with values computed independently of the library, read from the files in data/ (a test without its file is skipped)
add_executable(reference_tests
	reference/math_reference_tests.cpp
	reference/pdf_reference_tests.cpp
	reference/sidis_reference_tests.cpp
)

foreach(test_program unit_tests physics_tests reference_tests)
	# With the optimization flags of the dimuon program
	target_link_libraries(${test_program} PRIVATE dimuon::dimuon dimuon::optimization_flags GTest::gtest_main)
	# The support headers are included as "support/<file>.hpp"
	target_include_directories(${test_program} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
	target_compile_definitions(${test_program} PRIVATE DIMUON_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()

# Each test case is a CTest test, labeled by its program. The test cases are listed when the tests are run, not when they are built.
gtest_discover_tests(unit_tests PROPERTIES LABELS unit DISCOVERY_MODE PRE_TEST)
gtest_discover_tests(physics_tests PROPERTIES LABELS physics DISCOVERY_MODE PRE_TEST)
gtest_discover_tests(reference_tests PROPERTIES LABELS reference DISCOVERY_MODE PRE_TEST)
