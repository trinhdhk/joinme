# Declare association transformations with validation

`joinme_tf()` is a user-facing wrapper for transformation declarations
used by
[`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md) and
[`simulate_joinme()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme.md).
It validates channel names, normalises common shorthands, and returns a
structured object that can be passed directly as the `transforms`
argument.

Supported channels are `cv_total`, `cs_total`, `cv_mean`, `cs_mean`,
`cv_marker`, `cs_marker`, `corr`, and `vcov`. The `corr` channel refers
to off-diagonal entries of the subject-specific Cholesky-correlation
factor `K`, while `vcov` refers to those same off-diagonal `K` entries
together with the subject-specific standard deviations.

Supported shorthands per channel:

- `"identity"`, `"ispline"`, `"ispline_penalised"`,
  `"ispline_penalized"`, `"ispline_expit"`, `"ispline_expit_penalised"`,
  `"ispline_expit_penalized"`, `"ispline_exp_penalised"`,
  `"ispline_exp_penalized"`, `"pwlin"`, `"functional"`

- a formula such as `~ log1p(x)`, which is interpreted as
  `list(type = "functional", expr = ~ log1p(x))`

- a full named list specification such as
  `list(type = "ispline", knots = ..., coeff = ..., degree = 3)`

## Usage

``` r
joinme_tf(..., .validate = TRUE)

jm_tf(..., .validate = TRUE)
```

## Arguments

- ...:

  Named channel specifications.

- .validate:

  Logical; if `TRUE` (default), validate the resulting transform
  declarations immediately.

## Value

An object of class `joinme_tf`.

## Examples

``` r
tf <- joinme_tf(
  cv_total = "identity",
  corr = ~ -x,
  vcov = "identity",
  cv_marker = list(type = "pwlin", x = c(-1, 0, 1), y = c(0.2, 1, 0.2))
)
print(tf)
#> Transformation specification for Joint Mixed Effects model
#>    channel    specification
#>   cv_total         identity
#>       corr functional (~-x)
#>       vcov         identity
#>  cv_marker  pwlin, points=3
```
