# WAIC for JoiNMe models

WAIC for JoiNMe models

## Usage

``` r
# S3 method for class 'JoiNMeFit'
waic(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...)
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
  [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html).

## Value

A `waic` object.
