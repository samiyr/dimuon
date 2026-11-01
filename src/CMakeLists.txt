option(DIMUON_NATIVE "Optimize for the processor of this machine (-march=native)" ON)

# Propagate warning flags to all targets defined after this (`-Wnrvo` exists only for C++, not for the Fortran tools)
add_compile_options(-Wall -Wconversion -Wextra -Wpedantic -Wuninitialized -Wshadow $<$<COMPILE_LANGUAGE:CXX>:-Wnrvo>)

# Add the ability to check that library header files build independently via `make verify_headers`
set(CMAKE_VERIFY_INTERFACE_HEADER_SETS ON)
add_custom_target(verify_headers)
add_dependencies(verify_headers all_verify_interface_header_sets)

# Add the library
add_subdirectory(libdimuon)
# Add the argument parser of the executables
add_subdirectory(libargparser)
# Add the executable
option(BUILD_DRIVER "Build the main driver executable" ON)
if (BUILD_DRIVER)
	add_subdirectory(driver)
endif()
# Add the tools, such as the decay grid generator
option(BUILD_TOOLS "Build the tools" OFF)
if (BUILD_TOOLS)
	add_subdirectory(tools)
endif()
# Add the tests (testing itself is enabled in the top-level CMakeLists.txt)
option(BUILD_TESTS "Build the tests" OFF)
if (BUILD_TESTS)
	add_subdirectory(tests)
endif()

# Install the pool controller to <directory>/opt/dimuon/pool_controller
install(DIRECTORY pool_controller DESTINATION opt/dimuon)