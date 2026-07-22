# Summary of JoiNMe dynamic prediction

Summarises a `JoiNMeDynPred` object

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
posterior_summary(object, ...)

# S3 method for class 'JoiNMeDynPred'
summary(object, ...)
```

## Arguments

- object:

  A `JoiNMeDynPred` object produced by
  [predict](https://rdrr.io/r/stats/predict.html).

- ...:

  Unused.

## Value

A `summary_JoiNMeDynPred` object containing tabular summaries.

## Details

The summary includes:

- an overview table of row/subject counts by process,

- median survival time per subject (draw-wise crossing of `S(t)=0.5`,
  summarised with estimate, uncertainty, and diagnostics),

- predicted id-level random effects per subject (estimate, uncertainty,
  interval, and diagnostics),

- predicted marker-by-id random effects and covariance per subject when
  marker covariance depends on id,

- a compact diagnostics table counting potential convergence/ESS issues.
