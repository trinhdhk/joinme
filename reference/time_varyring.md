# Declare a time-varying simulation covariate generator

Helper used inside `covariate_formulas` in
[`simulate_joinme()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme.md)
to create subject-specific stepwise covariates over time. Elsewhere,
this function does not evaluate anything. Only a specification list is
returned.

## Usage

``` r
time_varyring(fun, n_step = NULL, steps = NULL, ...)
```

## Arguments

- fun:

  Random generator function (for example `rnorm`, `rbinom`).

- n_step:

  Integer vector of length 1 or 2.

  - length 1: each subject receives exactly `n_step` step periods.

  - length 2: each subject draws its number of periods uniformly from
    `min(n_step):max(n_step)`. Ignored when `steps` is supplied.

- steps:

  Optional vector of positive integer time points where the step changes
  occur exactly. Defaults to `NULL`. Values greater than
  subject-specific stop times are ignored.

- ...:

  Additional arguments forwarded to `fun`.

## Value

Specification fed into
[`simulate_joinme()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme.md).
