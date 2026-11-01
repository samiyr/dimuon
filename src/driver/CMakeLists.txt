# Create the dimuon executable
add_executable(dimuon main.cpp)
# Link libdimuon with the optimization flags, and the argument parser
target_link_libraries(dimuon PRIVATE dimuon::dimuon dimuon::optimization_flags argparser)
# Install to <directory>/bin
install(TARGETS dimuon)
