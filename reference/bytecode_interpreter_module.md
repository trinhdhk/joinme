# Portable stack-bytecode interpreter

A model-agnostic interpreter for scalar transformation programmes
encoded in reverse Polish notation. This file deliberately depends only
on base R and the stats package. It does not inspect JoiNMe formulae,
fitted objects or Stan data, which allows the file to be moved into a
small standalone package without changing its execution contract.

The authoritative instruction protocol is documented in
`inst/stan/include/etc/bytecode/README.md` and mirrored by the Stan
module under `inst/stan/include/etc/bytecode`.
