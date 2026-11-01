Copied unmodified from https://github.com/Expander/polylogarithm at commit 5ed2ca977fff734a2c5ce6e0ccceda9b8faf422f:
`src/cpp/Li2.cpp`, `Li2.hpp`, `Li3.cpp`, `Li3.hpp`, `horner.hpp` and `log.hpp` (which the two sources include), and `LICENSE` (MIT).
Not copied: the other functions of the library, and its C and Fortran implementations.

`Li2.cpp` and `Li3.cpp` are not compiled on their own: `dimuon/utility/math.hpp` includes them like headers, inside an unnamed namespace.
