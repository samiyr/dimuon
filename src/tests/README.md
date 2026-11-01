# Tests

The tests of `libdimuon`. They are built only when CMake is configured with `-DBUILD_TESTS=ON` (see installation instructions in the top-level `README`).

## Running

The tests are run with CTest in the build directory, either all of them or those with a label:

```
cd build
ctest
ctest -L physics
```