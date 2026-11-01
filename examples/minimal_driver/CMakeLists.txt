cmake_minimum_required(VERSION 3.31)

project(minimal_driver LANGUAGES CXX)

# Find the installed dimuon package (point CMake to it with e.g. `-DCMAKE_PREFIX_PATH=<installdir>`)
find_package(dimuon REQUIRED)

add_executable(minimal_driver main.cpp)
# Link the library and the recommended optimization flags
target_link_libraries(minimal_driver PRIVATE dimuon::dimuon dimuon::optimization_flags)
