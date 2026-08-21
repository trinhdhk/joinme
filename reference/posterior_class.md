# Posterior class membership for a latent-progress fit

Returns posterior class probabilities for the natural allocation domains
of a
[`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md)
model. Subject, correlation and covariance classes are reported by
subject; marker classes are reported by marker.

## Usage

``` r
posterior_class(object, ...)

# S3 method for class 'JoiNMeMixFit'
posterior_class(object, draws = NULL, seed = 1, digits = 3, summary = TRUE)
```

## Arguments

- object:

  A `JoiNMeMixFit` object.

- draws:

  Optional number of posterior draws.

- seed:

  Seed used for reproducible posterior subsetting.

- digits:

  Number of decimal places used for posterior summaries.

- summary:

  If `TRUE`, return the posterior mean (`Estimate`), posterior standard
  deviation (`Est.Error`), interval bounds, available sampling
  diagnostics, and the maximum-probability class. If `FALSE`, return
  draw matrices.

## Value

A named list with `subject` and/or `marker` entries.
