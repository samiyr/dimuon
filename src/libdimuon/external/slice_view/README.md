Copied from https://github.com/bemanproject/sandbox-slice_view at commit e3e6130cded2765a2c30d2a96f754e6120268c46:
`include/beman/slice_view/slice_view.hpp` and `detail/concepts.hpp` (which the former includes), and `LICENSE` (Apache-2.0 WITH LLVM-exception).

One change: `slice_view.hpp` includes `"detail/concepts.hpp"` instead of `<beman/slice_view/detail/concepts.hpp>`, since the files are not below a `beman/slice_view` directory here.
