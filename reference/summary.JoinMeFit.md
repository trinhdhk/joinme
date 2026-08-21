# Summary of a JoiNMe object

Builds posterior summaries for the longitudinal process, survival
process, association terms, and optional covariance blocks. For
survival, the summary reports `survival_process` when the event model
contains non-intercept covariates beyond association features.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
summary(
  object,
  draws = NULL,
  seed = .Random.seed[[1]],
  digits = 3,
  include_corr = TRUE,
  ...
)

# S3 method for class 'JoiNMeFit'
posterior_summary(
  object,
  what = c("model", "fixef", "ranef", "coef"),
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  A JoiNMe fit object.

- draws:

  Number of draws to use for summaries.

- seed:

  Random seed for subsetting draws.

- digits:

  Number of digits to round summary values.

- include_corr:

  Logical; include covariance summaries.

- ...:

  Unused.

- what:

  Posterior view to return. `"model"` gives the complete fitted model
  summary. `"fixef"`, `"ranef"`, and `"coef"` provide the common
  reporting layer used by the corresponding high-level methods.

- summary:

  Logical. For a coefficient view, `TRUE` returns posterior summaries
  and `FALSE` returns its structured draw-level extraction.

## Value

A `summary_JoiNMeFit` object. Its `tables` element includes posterior
association coefficients, a compact marker-weight table with common
locations and within-set spreads, transform parameters, and, when an
ordered piecewise-linear association is active, `piecewise_ordinates`
containing the relative log-hazard and hazard-ratio contribution at
every knot. Declared marker-weight offsets are retained in
`metadata$marker_weight_offsets` and printed above the posterior tables.

## Details

For coefficient views, this is the reporting layer between
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
and the conventional
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html) methods. The high-level
methods delegate here; requests with `summary = FALSE` continue to the
component-aware draw representation owned by
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md).
