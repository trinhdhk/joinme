# Transform Specification Helper

Transform can be specified in four modes:

- Mode 0: Identity (no transformation)

- Mode 1: Functional (arbitrary nested functions via parser)

- Mode 2: I-spline basis (monotonic spline)

- Mode 2 (penalised): I-spline coefficients estimated with smoothness
  penalty

- Mode 4: I-spline basis evaluated on `expit(x)` for a bounded,
  numerically stable spline input domain

- Mode 3: Increasing ordered piecewise-linear association

- Mode 7: Decreasing ordered piecewise-linear association

Each transformation term (CV_total, CS_total, CV_mean, CS_mean,
CV_marker, CS_marker, corr, vcov) can have its own independent
specification, enabling flexibility in model building. For
covariance-style channels, `corr` transformations act on off-diagonal
entries of the subject-specific Cholesky-correlation factor `K`, while
`vcov` transformations act on those same off-diagonal `K` entries
together with the subject-specific standard deviations.

## Usage

``` r
build_standata_transforms(
  transform_list = NULL,
  default_mode = 0,
  n_corr_components = NULL,
  n_vcov_components = NULL
)
```

## Arguments

- transform_list:

  List with elements cv_total, cs_total, cv_mean, cs_mean, cv_marker,
  cs_marker, corr, vcov, each specifying a transformation. See details.

- default_mode:

  Default transformation mode if not specified.

- n_corr_components:

  Optional non-negative number of correlation components. It determines
  the row count of component-specific coefficient matrices.

- n_vcov_components:

  Optional non-negative number of covariance components. It determines
  the row count of component-specific coefficient matrices.

## Value

List of standata entries for transformation parameters:

- tf_mode\_\*

- functional_ops\_\* and const_data\_\* (if mode 1)

- estimate_iota_intercept\_\* and estimate_iota_slope\_\* (fit-only
  functional affine-shift counts)

- functional_iota_intercept_idx\_\* and functional_iota_slope_idx\_\*
  (per-op affine-shift indices)

- knots\_\* and coeff\_\* and spline_degree\_\* (if mode 2, 3, or 4)

## Details

Provides utilities to build and manage transformation specifications for
the JoiNMe Stan model's association term.

Each element of transform_list should be a list with:

- type: "identity", "functional", "ispline", "ispline_penalised" (or
  alias "ispline_penalized"), "ispline_expit", "ispline_expit_penalised"
  (aliases `"ispline_exp_penalised"`, `"ispline_exp_penalized"`,
  `"ispline_expit_penalized"`), or "pwlin"

- Additional fields depend on type:

  - functional: expr (quosure, formula, quoted expression, or string)

  - ispline: knots (vector), coeff (vector), degree (int)

  - ispline_penalised: knots or n_knots, degree, lambda, optional x,
    optional y, optional weights

    - if y is supplied, JoiNMe fits the monotone spline to the training
      pairs (x, y) in R using the same anchored endpoint convention as
      the Stan-estimated path (first coefficient = 0, last coefficient =
      1)

    - if y is omitted, Stan estimates the monotone spline coefficients
      directly and lambda controls smoothness

  - ispline_expit / ispline_expit_penalised: identical to the I-spline
    variants above, except the raw association feature is first mapped
    through `plogis(x)` and the spline basis is then evaluated on that
    bounded expit-scale input. User-supplied training `x` values and
    explicit `knots` for these transform types must therefore already be
    specified on the expit scale in `[0, 1]`.

  - pwlin: knots (or cutpoints/x) and direction; legacy y is accepted
    but does not determine the fitted ordinates

Defaults and minimal examples:

- identity (default if term omitted): `list(type = "identity")`

- functional: `list(type = "functional", expr = ~ log1p(x))`

- ispline:
  `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3))`
  with `degree = 3` by default when omitted.

- ispline_penalised (plug-in fit):
  `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`

- ispline_penalised (Stan-estimated):
  `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`

- ispline_expit:
  `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1))`

- ispline_expit_penalised (plug-in fit):
  `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0, 1, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`

- ispline_expit_penalised (Stan-estimated):
  `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`

- pwlin:
  `list(type = "pwlin", knots = c(-2, -1, 0, 1, 2), direction = "increasing")`

Default values used internally:

- omitted term -\> identity mode (`default_mode = 0`),

- `ispline`: `degree = 3` if omitted,

- `ispline_penalised`: `n_knots = 6`, `degree = 3`, `lambda = 1` if
  omitted,

- `pwlin`: `direction = "increasing"` if omitted,

- `functional`: no additional defaults beyond its required fields.

For `ispline_expit*` declarations, the internal spline basis lives on
the bounded interval \$(0, 1)\$ after applying
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) to the raw
association feature. User-supplied training `x` values and explicit
`knots` should be on that same bounded scale so the fitted spline domain
matches the runtime basis. This often yields better-conditioned knot
placement and smoother optimisation for steep nonlinear transforms.

For user-facing declarations, prefer `joinme_tf(...)`, which validates
and normalises channel specifications before they are passed here.

## Usage

Use `build_standata_transforms()` to construct the data list entries,
and `validate_transforms()` to verify consistency before sampling.
