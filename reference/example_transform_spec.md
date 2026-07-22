# Example: Build Transform Specification from Formula String

Dummy function to print help this help page. See
[`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md)
for details on how to construct a transform for the association term.

## Usage

``` r
example_transform_spec()
```

## Examples

``` r
# Example 1: Functional transformation (bytecode-backed)
spec_cv <- list(
  type = "functional",
  expr = ~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
)

# Example 2: Piecewise-linear
spec_cs <- list(
  type = "pwlin",
  x = c(-2, -1, 0, 1, 2),
  y = c(0.1, 0.3, 1.0, 0.3, 0.1)  # smooth bump
)

# Example 3: I-spline (monotonic)
# - knots define basis locations (first/last are boundary knots)
# - coeff defines the monotone shape directly
spec_corr <- list(
  type = "ispline",
  knots = c(-1, 0, 1),
  coeff = c(0, 0.5, 1, 1.2, 1.5),
  degree = 3
)

# Example 4: Penalised monotone I-spline (fit in R)
# - x is the input scale of the raw association feature
# - y is the desired transformed output at each x
# - lambda controls smoothness (higher = smoother)
spec_corr_pen <- list(
  type = "ispline_penalised",
  x = seq(-2, 2, length.out = 50),
  y = exp(seq(-2, 2, length.out = 50)),
  n_knots = 6,
  degree = 3,
  lambda = 1.0,
  direction = "increasing"
)

# Example 5: Penalised monotone I-spline (Stan-estimated coefficients)
# - omit y so Stan learns the monotone shape directly
# - if knots are omitted, n_knots = 6 and x quantiles define them
spec_corr_pen_stan <- list(
  type = "ispline_penalised",
  x = seq(-2, 2, length.out = 50),
  n_knots = 6,
  degree = 3,
  lambda = 1.0
)

# Example 6: Penalised monotone I-spline on expit(x)
spec_corr_expit <- list(
  type = "ispline_expit_penalised",
  x = seq(0.02, 0.98, length.out = 50),
  y = seq(0.02, 0.98, length.out = 50)^0.75,
  n_knots = 6,
  degree = 3,
  lambda = 1.0,
  direction = "increasing"
)

# Combine and build standata
transforms <- joinme_tf(
  cv_total = spec_cv,
  cs_total = spec_cs,
  corr = spec_corr
)
standata_tf <- build_standata_transforms(transforms)
#> Error in build_standata_transforms(transforms): could not find function "build_standata_transforms"
```
