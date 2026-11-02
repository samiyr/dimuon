This code is for computing dimuon production cross sections in neutrino-nucleus collisions in the SIDIS framework, based on the following papers:
- [JHEP 09 (2024) 043](https://doi.org/10.1007/JHEP09(2024)043)
- [JHEP 01 (2026) 112](https://doi.org/10.1007/JHEP01(2026)112)
- [JHEP 06 (2026) 254](https://doi.org/10.1007/JHEP06(2026)254).

# Components

The code consists of several components, which are listed here. Each component has more documentation within it.

- `libdimuon`: header-only library that actually does the computations
- `driver`: driver executable, which calls into `libdimuon`, that can be run to produce cross section or interpolation table values
- `pool_controller`: a Python controller for running `dimuon` on a cluster
- `libargparser`: supporting header-only library for parsing command-line arguments
- `tools/decay_grid_generator`: tool for calculating decay grids
- `tools/bkk05_lhapdf_grid_generator`: tool for converting the [`bkk05`](https://doi.org/10.1103/PhysRevD.74.037502) fragmentation function grids into the LHAPDF format
- `tools/kkks08_lhapdf_grid_generator`: tool for converting the [`kkks08`](https://doi.org/10.1016/j.nuclphysb.2008.02.015) fragmentation function grids into the LHAPDF format
- `tools/interpolation_tables`: Python tools for reading interpolation tables generated from `dimuon` and converting them to the [PineAPPL](https://zenodo.org/badge/latestdoi/248306479) format
- `tests`: various tests of `libdimuon`

# Installation

The dimuon code is installed with `cmake`, which does four things by default:
1. Builds the `dimuon` executable and copies it to `<installdir>/bin`
2. Installs the library header files into `<installdir>/include`
3. Installs the decay grid files into `<installdir>/share/dimuon/decaygrids`
4. Installs the example driver (without building it) into `<installdir>/share/dimuon/examples`

Since `libdimuon` is a header-only library, the header files are simply copied and not built. The library can only be built when it's included in another target (a driver).

The library code is written in `c++23` and needs a modern compiler. The code works on at least `GCC 15.2` and newer.

## Dependencies

Before installing, the following dependencies have to be installed:
- [`Cuba`](https://feynarts.de/cuba/)
- [`LHAPDF`](https://www.lhapdf.org)
- [`GSL`](https://www.gnu.org/software/gsl/)

For CMake to find these, it might be necessary to point the `CUBA_ROOT`, `LHAPDF_ROOT`, and `GSL_ROOT` variables to the respective installations so that CMake can find these libraries (see [`find_package`](https://cmake.org/cmake/help/latest/command/find_package.html) documentation for more).

After installing these libraries, you need to set the proper environment variables. Regardless of the build system, `LD_LIBRARY_PATH` has to contain the appropriate paths for the dynamic linker to find the libraries at runtime. If you get `shared library not found` errors related to these libraries when running an executable, check `LD_LIBRARY_PATH`.

## Installing the library

To install, first run the configuration command in this directory
```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<installdir>
```
where `<installdir>` is the location where the library will be installed. If necessary, pass on flags like `-DCUBA_ROOT=<cubadir>`. Other flags that can be used with the above command:
- `-DDIMUON_NATIVE=OFF`: removes the `-march=native` flag (on by default)
- `-DBUILD_DRIVER=OFF`: skips building the `dimuon` executable, in which case there's nothing to build in the next step (on by default)
- `-DBUILD_TOOLS=ON`: builds the executables in the `tools` directory (off by default)
- `-DBUILD_TESTS=ON`: builds the test executables (off by default)

After configuring, run the build command to build the `dimuon` executable (and other executables if enabled)
```
cmake --build build
```
Use the `-j` flag to build in parallel.

Finally, run the install command to copy all files to `<installdir>`:
```
cmake --install build
```

## Optimization flags

The library cannot have any optimization flags enabled by default, and instead have to be set on the target that includes `libdimuon`. Recommended flags are exported as `dimuon::optimization_flags` for easy access via CMake.

## Drivers

Again, the library itself is header-only and cannot be built on its own, so there will be no `lib` directory in `<installdir>` and no `libdimuon.so` file. The library will be built only when it's included in a driver executable.

A minimal example of a driver can be found in `examples/minimal_driver`. The `dimuon` executable is an example of a more complete driver that can be used for obtaining actual data. The built `dimuon` executable can be found in `<installdir>/bin/dimuon`.

To verify that the library builds properly without including it in a driver, the command
```
cmake --build build --target verify_headers
```
can be run. This command compiles each translation unit (= header file) independently.

### Compiling a driver directly

As usual, you need to point the compiler to the `<installdir>/include` and `<installdir/include/dimuon/external>` directory (likely with `-I`, depending on your compiler), and to any other necessary paths for the dependencies. As an example, the minimal driver can be built directly using `gcc` with the command
```
cd <installdir>
g++ share/dimuon/examples/minimal_driver/main.cpp -std=c++23 \
-O3 -DNDEBUG -fno-math-errno -fno-trapping-math -fassociative-math -fno-signed-zeros -march=native \
-Iinclude -Iinclude/dimuon/external \
-I<path_to_lhapdf>/include -I<path_to_cuba>/include \
-L<path_to_lhapdf>/lib -L<path_to_cuba>/lib \
-lLHAPDF -lcuba -lgsl -lgslcblas \
-o share/dimuon/examples/minimal_driver/minimal_driver.out
```
The second line contains the optimization flags: `-O3 -DNDEBUG` are the flags CMake uses for a `Release` build, and the rest are the recommended flags exported as `dimuon::optimization_flags`. The includes of the external dependencies on the following lines will depend on how they are configured in the environment, and whether they are already included in e.g. `CPATH`. In the example above, the GSL include and linker paths are missing because they are already included through `CPATH`, `LIBRARY_PATH`, and `LD_LIBRARY_PATH`.

### Compiling a driver with CMake

If necessary, point the `CUBA_ROOT`, `LHAPDF_ROOT`, and `GSL_ROOT` variables to the respective installations so that CMake can find these libraries (see [`find_package`](https://cmake.org/cmake/help/latest/command/find_package.html) documentation for more).

As an example, the minimal driver also contains a `CMakeLists.txt`. It finds the installed library and links against it with the appropriate optimization flags. It can be built with the commands
```
cd <installdir>/share/dimuon/examples/minimal_driver
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=<installdir>
cmake --build build
```
which creates the executable `build/minimal_driver`. Note that `-DCMAKE_BUILD_TYPE=Release` is necessary for the optimizations to be enabled.

# Available data

Some of the data files are publicly available via [this link](https://a3s.fi/swift/v1/jyy2580-dimuon-pub/). The simplest way to use that link is to open it in a browser, which will give you a list of file paths. The link to an actual file is then obtained by appending a file path to the link above.

The most relevant files are the ones in the `dimuon-data/pineappl` directory, which are the dimuon production interpolation tables in the [PineAPPL](https://zenodo.org/badge/latestdoi/248306479) format. For more information on how to use these PineAPPL grids, refer to the `README.md` in `tools/interpolation_tables`. 