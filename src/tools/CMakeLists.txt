# Add the decay grid generator
option(SKIP_BUILD_DECAY_GRID_GENERATOR "Skip the build and install of the decay grid generator" OFF)
if (NOT SKIP_BUILD_DECAY_GRID_GENERATOR)
	add_subdirectory(decay_grid_generator)
endif()

# Add the LHAPDF grid generator of the KKKS08 fragmentation functions
option(SKIP_BUILD_KKKS08_LHAPDF_GRID_GENERATOR "Skip the build and install of the LHAPDF grid generator of the KKKS08 fragmentation functions" OFF)
if (NOT SKIP_BUILD_KKKS08_LHAPDF_GRID_GENERATOR)
	enable_language(Fortran)
	add_subdirectory(kkks08_lhapdf_grid_generator)
endif()

# Add the LHAPDF grid generator of the BKK05 fragmentation functions
option(SKIP_BUILD_BKK05_LHAPDF_GRID_GENERATOR "Skip the build and install of the LHAPDF grid generator of the BKK05 fragmentation functions" OFF)
if (NOT SKIP_BUILD_BKK05_LHAPDF_GRID_GENERATOR)
	enable_language(Fortran)
	add_subdirectory(bkk05_lhapdf_grid_generator)
endif()

# Install the gc coefficient generator to <directory>/opt/dimuon/gc_coefficient_generator
option(SKIP_INSTALL_GC_GENERATOR "Skip the installation of the gc SIDIS coefficient generator" OFF)
if (NOT SKIP_INSTALL_GC_GENERATOR)
	install(DIRECTORY gc_coefficient_generator DESTINATION opt/dimuon)
endif()

# Install the interpolation table scripts to <directory>/opt/dimuon/interpolation_tables
option(SKIP_INSTALL_INTERPOLATION_TABLES "Skip the installation of the interpolation table scripts" OFF)
if (NOT SKIP_INSTALL_INTERPOLATION_TABLES)
	install(DIRECTORY interpolation_tables DESTINATION opt/dimuon)
endif()