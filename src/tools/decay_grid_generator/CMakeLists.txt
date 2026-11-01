# The grids are computed in parallel with OpenMP
find_package(OpenMP REQUIRED)

# Create the decay grid generator executable
add_executable(decay_grid_generator main.cpp)
# Link libdimuon with the optimization flags, the argument parser, and OpenMP
target_link_libraries(decay_grid_generator PRIVATE dimuon::dimuon dimuon::optimization_flags argparser OpenMP::OpenMP_CXX)
# Install to <directory>/bin
install(TARGETS decay_grid_generator)
