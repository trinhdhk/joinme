# Log-likelihood summary

Log-likelihood summary

## Usage

``` r
# S3 method for class 'JoiNMeFit'
log_lik(object, what = c("long", "surv", "total"), draws = NULL, seed = 1, ...)
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

  Unused.

## Value

A matrix

## See also

[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
