# Extract predicted random effects from dynamic predictions

Returns predicted random effects from a `JoiNMeDynPred` object.
Marker-by-id random effects are available only when marker covariance is
configured to be subject-dependent.

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
ranef(object, ...)
```

## Arguments

- object:

  A `JoiNMeDynPred` object.

- ...:

  Unused.

## Value

A named list containing random-effects summary tables.
