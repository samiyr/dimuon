# Create the BKK05 LHAPDF grid generator executable from the generator and the original code (bkk05_D3.f)
add_executable(bkk05_lhapdf_grid_generator bkk05_D3_grid_generator.f90 bkk05_D3.f)
# Compile the original code
set_property(SOURCE bkk05_D3.f PROPERTY COMPILE_OPTIONS -Wno-pedantic -std=legacy -w)
# Install to <directory>/bin
install(TARGETS bkk05_lhapdf_grid_generator)
