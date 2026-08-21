# Functional-transformation submodel

This directory describes transformations applied to a longitudinal association
feature before its hazard coefficient is used.

A transformation may be a fixed bytecode expression, a monotone spline, an
ordered piecewise-linear curve, or an expression with fitted affine intercept
and slope. Unit-span shape constraints distinguish transformation shape from
association scale.

The fitted affine coefficients are named `iota_intercept_*` and
`iota_slope_*`. They modify a feature inside `tf()` before marker averaging;
they are distinct from the `alpha_*` association slope that multiplies the
completed feature in the event log hazard.

## Read the blocks in order

1. [data](data/README.md) — modes, bytecode, constants, knots and fitted
   affine-term flags.
2. [transformed data](transformed_data/README.md) — validated indexing and
   parameter-free preparation.
3. [parameters](parameters/README.md) — shape, simplex and affine unknowns.
4. [transformed parameters](transformed_parameters/README.md) — unit-span
   curves and effective affine terms.
5. [model](model/README.md) — shape, smoothness and affine prior contributions.
6. [generated quantities](generated_quantities/README.md) — fitted affine
   outputs and trajectory ownership.

[Back to the master submodel map](../README.md).
