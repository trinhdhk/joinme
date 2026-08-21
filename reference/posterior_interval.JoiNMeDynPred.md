# Central posterior credible intervals for dynamic predictions

Flattens the stored dynamic-prediction draws using
[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md)
and computes central credible intervals with
[`rstantools::posterior_interval()`](https://mc-stan.org/rstantools/reference/posterior_interval.html).
Composite parameter names retain the subject, marker, prediction scale,
and evaluation-time identities.

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
posterior_interval(
  object,
  prob = 0.9,
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  ...
)
```

## Arguments

- object:

  A `JoiNMeDynPred` dynamic-prediction object. Latent-class predictions
  inherit this method through `JoiNMeMixDynPred`.

- prob:

  A single number strictly between zero and one giving the posterior
  probability contained in the interval. The default is `0.9`, following
  `rstantools`.

- variables:

  Optional character vector selecting friendly parameter names when
  `what = "model"`.

- regex:

  Logical; when `TRUE`, interpret `variables` as regular expressions.
  This argument applies only to `what = "model"`.

- draws:

  Optional number of posterior draws to retain before computing the
  intervals.

- seed:

  Integer seed used when posterior draws are subsampled.

- ...:

  Additional arguments passed to
  [`rstantools::posterior_interval()`](https://mc-stan.org/rstantools/reference/posterior_interval.html).

## Value

A numeric matrix with one row per selected prediction quantity and two
probability-labelled interval columns.
