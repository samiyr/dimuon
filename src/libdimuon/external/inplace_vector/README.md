Copied from https://github.com/bemanproject/inplace_vector at commit c7fe76da307f0d6fda5023b2d2ebef7e3f7b18d1:
`include/beman/inplace_vector/inplace_vector.hpp` and `config.hpp` (which the former includes), and `LICENSE` (Apache-2.0 WITH LLVM-exception).

One change: `inplace_vector.hpp` includes `"config.hpp"` instead of `<beman/inplace_vector/config.hpp>`, since the files are not below a `beman/inplace_vector` directory here.
