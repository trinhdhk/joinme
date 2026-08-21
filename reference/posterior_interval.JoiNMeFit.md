# Central posterior credible intervals for JoiNMe models

Computes central posterior credible intervals through
[`rstantools::posterior_interval()`](https://mc-stan.org/rstantools/reference/posterior_interval.html).
The method uses the same scientific coefficient views as
[`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html):
the complete fitted posterior, population-level coefficients,
group-specific deviations, or their combined coefficients.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
posterior_interval(
  object,
  prob = 0.9,
  what = c("model", "fixef", "ranef", "coef"),
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  ...
)
```

## Arguments

- object:

  A fitted `JoiNMeFit` object. Latent-class fits inherit this method
  through `JoiNMeMixFit`.

- prob:

  A single number strictly between zero and one giving the posterior
  probability contained in the interval. The default is `0.9`, following
  `rstantools`.

- what:

  Posterior view to interval: `"model"`, `"fixef"`, `"ranef"`, or
  `"coef"`.

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

For `what = "model"` or `what = "fixef"`, a numeric matrix with one row
per term and two probability-labelled columns. For `what = "ranef"` or
`what = "coef"`, a nested list of data frames retaining the coefficient
identifiers and containing the same two probability-labelled columns.

## Details

For `what = "model"`, the returned matrix has one row per friendly
posterior parameter name. The `variables` argument selects those names
after JoiNMe has translated the Stan coordinates into their reported
statistical terms.

The `fixef` view is also rectangular and therefore returns an interval
matrix. The `ranef` and `coef` views contain several statistically
distinct coefficient tables. Their nested list structure is retained,
while each draw-level `value` column is replaced by the two interval
limits. This keeps subject, marker, event, distributional-family,
association-term, and covariance identities explicit.

Every interval is calculated by the default matrix method of
[`rstantools::posterior_interval()`](https://mc-stan.org/rstantools/reference/posterior_interval.html).
Consequently, `prob` is the total posterior probability contained
between the two central quantiles.
