# Build named conditional-effects profiles

`make_conditions()` is a thin wrapper around
[`brms::make_conditions()`](https://paulbuerkner.com/brms/reference/make_conditions.html)
so the resulting labelled condition tables can be passed directly to
[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md).
Conditioning belongs to the conditional- effects estimand; neither
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md)
nor
[`plot.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/plot.JoinMeFit.md)
accepts a `condition` argument.

## Usage

``` r
make_conditions(x, ...)
```

## Arguments

- x:

  A data frame containing the baseline covariates used to define the
  conditioning rows.

- ...:

  Additional arguments passed to
  [`brms::make_conditions()`](https://paulbuerkner.com/brms/reference/make_conditions.html).

## Value

A data frame with one row per requested condition and a `cond__` column
containing the profile labels used in conditional-effects facets.

## See also

[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md),
[`brms::make_conditions()`](https://paulbuerkner.com/brms/reference/make_conditions.html)
