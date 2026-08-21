# Predict from a latent-progress mixture

Runs the established dynamic prediction calculation while carrying the
fitted component probability, location and scale into the priors for
newly sampled latent effects. The new subject's history therefore
informs a posterior probability vector over the fitted classes. Combined
compatible blocks retain one shared allocation. With
`reuse_fitted_re = TRUE`, the fitted subject's realised effects and
fitted allocation probabilities are retained draw by draw instead of
evaluating a new-subject allocation.

## Usage

``` r
# S3 method for class 'JoiNMeMixFit'
predict(
  object,
  newdataLong,
  newdataEvent = NULL,
  process = c("longitudinal", "event"),
  reuse_fitted_re = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- newdataLong:

  Data frame containing longitudinal histories for one or more subjects.
  Must contain columns for id, time, marker, and response variables as
  specified in the original model formula.

- newdataEvent:

  Data frame containing event information and covariates used in the
  survival model. It may contain one row per subject or interval-split
  rows (`Surv(start, stop, status)` layout). Left- and interval-censored
  survival encodings (`type = "left"`, `type = "interval2"`) are
  accepted. When interval rows are supplied, dynamic prediction uses the
  latest row per subject for event-side covariates. This argument may be
  omitted for a longitudinal-only fit; neutral event rows are then
  constructed internally.

- process:

  Character vector specifying which predictions to compute. Options:
  "longitudinal" (future trajectory), "event" (conditional survival
  probability). Default: both.

- reuse_fitted_re:

  Logical. If `TRUE`, every identifier in the supplied data must have
  occurred during fitting. Prediction reuses that subject's paired
  posterior random effects and covariance draws, without estimating a
  second set of random effects. The default, `FALSE`, conditions newly
  sampled effects on the supplied longitudinal history for dynamic
  prediction.

- ...:

  Additional arguments passed to
  [predict](https://rdrr.io/r/stats/predict.html).

## Value

A `JoiNMeMixDynPred` object inheriting from `JoiNMeDynPred`. Conditional
allocation draws are retained in `draws$posterior_class` and may be
summarised with
[`posterior_class()`](https://trinhdhk.github.io/joinme/reference/posterior_class.md).
