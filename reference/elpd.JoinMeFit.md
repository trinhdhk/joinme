# ELPD summary for JoiNMe models

ELPD summary for JoiNMe models

## Usage

``` r
# S3 method for class 'JoiNMeFit'
elpd(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- what:

  Character; which component to return: `"long"`, `"surv"`, or
  `"total"`.

- draws:

  Optional number of posterior draws to subset.

- seed:

  Random seed for draw subsetting.

- ...:

  Additional arguments passed to
  [`loo::elpd()`](https://mc-stan.org/loo/reference/elpd.html).

## Value

A data frame with ELPD and standard error.
