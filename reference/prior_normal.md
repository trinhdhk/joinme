# Declare a Normal coefficient prior

Creates a JoiNMe prior declaration without evaluating a density in R.
The location and scale may be scalars or vectors; scalar values are
recycled to the number of coefficients in the selected model block.

The `prior_` prefix is intentional. It avoids masking the short
constructor names exported by other Bayesian modelling packages.

## Usage

``` r
prior_normal(mu = 0, scale = 1)
```

## Arguments

- mu:

  Prior location or vector of coefficient-specific locations.

- scale:

  Positive prior scale or vector of coefficient-specific scales. Numeric
  lists, such as `list(1, 1, 2)`, are accepted and flattened in their
  supplied order. The same convention applies to `mu`.

## Value

A `joinme_prior_spec` object for use inside
[`jm_prior()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
