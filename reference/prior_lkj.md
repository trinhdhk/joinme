# Declare an LKJ correlation prior

Declare an LKJ correlation prior

## Usage

``` r
prior_lkj(eta = 1)
```

## Arguments

- eta:

  Positive LKJ concentration. `eta = 1` is uniform over correlation
  matrices; values above one favour correlations nearer zero.

## Value

A `joinme_lkj_prior` object for use inside
[`jm_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
