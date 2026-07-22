# Summary of a JoiNMe object

Builds posterior summaries for the longitudinal process, survival
process, association terms, and optional covariance blocks. For
survival, the summary reports `survival_process` when the event model
contains non-intercept covariates beyond association features.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
summary(object, draws = NULL, seed = 1, digits = 3, include_corr = TRUE, ...)

# S3 method for class 'JoiNMeFit'
posterior_summary(object, ...)
```

## Arguments

- object:

  A JoiNMe fit object.

- draws:

  Number of draws to use for summaries.

- seed:

  Random seed for subsetting draws.

- digits:

  Number of digits to round summary values.

- include_corr:

  Logical; include covariance summaries.

- ...:

  Unused.

## Value

A `summary_JoiNMeFit` object
