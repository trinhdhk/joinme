# Declare a Student-t coefficient prior

Declare a Student-t coefficient prior

## Usage

``` r
prior_student_t(df = 6, mu = 0, scale = 1)
```

## Arguments

- df:

  Positive fixed degrees of freedom. A single value is shared by the
  complete parameter block. Marker-weight departures instead use the
  bare family name `"student_t"`, whose degrees of freedom are fitted.

- mu:

  Prior location or vector of coefficient-specific locations.

- scale:

  Positive prior scale or vector of coefficient-specific scales. Numeric
  lists, such as `list(1, 1, 2)`, are accepted and flattened in their
  supplied order. The same convention applies to `mu`.

## Value

A `joinme_prior_spec` object for use inside
[`jm_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
