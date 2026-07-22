# Extract predicted covariance summaries

Returns per-subject covariance summaries for marker-by-id random effects
from a `JoiNMeDynPred` object. This method is available only when marker
covariance depends on id.

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
vcov(object, ...)
```

## Arguments

- object:

  A `JoiNMeDynPred` object.

- ...:

  Unused.

## Value

A named list containing covariance summary tables.
