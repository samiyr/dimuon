# Create the KKKS08 LHAPDF grid generator executable from the generator and the original code (kkks08_D.f, fint.f)
add_executable(kkks08_lhapdf_grid_generator kkks08_grid_generator.f90 kkks08_D.f fint.f)
# Compile the original code, which is not modified here, without warnings
set_property(SOURCE kkks08_D.f fint.f PROPERTY COMPILE_OPTIONS -w)
# Install to <directory>/bin
install(TARGETS kkks08_lhapdf_grid_generator)
