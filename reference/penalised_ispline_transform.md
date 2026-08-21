# Penalised Monotone I-spline Transform Spec

Fits a monotone I-spline transformation with a smoothness penalty and
returns a transform specification compatible with
[`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md).

This helper is the plug-in constructor for `type = "ispline_penalised"`.
Its defaults are `n_knots = 6`, `degree = 3`, `lambda = 1.0`,
`weights = NULL` (equal weights), and `diff_order = 2`. Returned
coefficients follow the same anchored convention as the Stan- estimated
path: increasing splines run from `0` to `1`, decreasing splines run
from `1` to `0`, and interior coefficients stay monotone in the chosen
direction. Example:
`penalised_ispline_transform(x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)))`

## Usage

``` r
penalised_ispline_transform(
  x,
  y,
  knots = NULL,
  n_knots = 6,
  degree = 3,
  lambda = 1,
  direction = NULL,
  weights = NULL,
  diff_order = 2
)

penalized_ispline_transform(...)
```

## Arguments

- x:

  Numeric vector of input values on the raw feature scale to transform.

- y:

  Numeric vector of target transformed values at `x` (same length as
  `x`).

- knots:

  Optional numeric vector of knots (including boundary knots).

- n_knots:

  Integer. If knots are not provided, number of knots to use (including
  boundary knots). Default 6.

- degree:

  Integer spline degree (default 3).

- lambda:

  Non-negative smoothness penalty weight (default 1.0). Larger values
  produce smoother fitted transforms; smaller values allow more local
  curvature.

- direction:

  Monotone direction. Accepted values are `"increasing"` and
  `"decreasing"`. When omitted, the direction is inferred from the
  supplied `(x, y)` pairs.

- weights:

  Optional non-negative weights (same length as x).

- diff_order:

  Integer difference order for penalty (default 2).

- ...:

  Arguments passed to `penalised_ispline_transform()`.

## Value

A list suitable for
[`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md).
