# Posterior predictive check for longitudinal outcomes

Posterior predictive check for longitudinal outcomes

## Usage

``` r
# S3 method for class 'JoiNMeFit'
pp_check(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  ci_level = 0.95,
  n_samples = 200,
  seed = 123,
  plot = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- newdataLong:

  New longitudinal data frame for predictions.

- newdataEvent:

  New event data frame for predictions.

- ci_level:

  Numeric; credible interval level (default 0.95).

- n_samples:

  Integer; number of posterior samples to use (default 200).

- ...:

  Additional arguments.

## Value

A list with observation-level summaries and overall diagnostics.

## Details

This crap is still under development.
