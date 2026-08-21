# Conditional class membership after dynamic prediction

Summarises the class probabilities obtained after conditioning a new
subject's random effects on their observed longitudinal history. For a
combined subject and covariance-regression mixture, the probabilities
refer to their one shared allocation.

## Usage

``` r
# S3 method for class 'JoiNMeMixDynPred'
posterior_class(object, draws = NULL, seed = 1, digits = 3, summary = TRUE)
```

## Arguments

- object:

  A `JoiNMeMixDynPred` object returned by
  [`predict.JoiNMeMixFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeMixFit.md).

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

A `JoiNMePosteriorClass` list. With `summary = TRUE`, the `subject` and
`marker` entries are posterior summary tables. With `summary = FALSE`,
the subject-indexed matrices and arrays are returned.
