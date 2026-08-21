# Declare a Laplace coefficient prior

The Laplace distribution is Stan's double-exponential distribution. Its
`scale` is the exponential-decay scale, not its standard deviation; the
standard deviation is `sqrt(2) * scale`.

## Usage

``` r
prior_laplace(mu = 0, scale = 1)
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
