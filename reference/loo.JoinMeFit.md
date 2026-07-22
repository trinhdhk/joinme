# LOO-CV for JoiNMe models

LOO-CV for JoiNMe models

## Usage

``` r
# S3 method for class 'JoiNMeFit'
loo(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...)
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
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object.
