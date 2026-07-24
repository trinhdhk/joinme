# Time-dependent area under the ROC curve

Estimate the cumulative/dynamic area under the receiver operating
characteristic curve for the survival process at one or more landmark
times. Cases experience the requested cause after the landmark and by
the horizon; controls remain event-free beyond the horizon. Subjects
censored before the horizon are excluded because their case/control
state is unknown.

Posterior mean dynamic risks are obtained through
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md),
so the calculation automatically uses the complete fitted hazard:
baseline hazard, event covariates, marker weights, and any identity,
functional, monotone-spline, or ordered piecewise-linear association
transform.

## Usage

``` r
auc(object, ...)

# S3 method for class 'JoiNMeFit'
auc(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  time_start = 0,
  time_horizon = NULL,
  Dt = NULL,
  cause = 1,
  n_samples = 200,
  seed = 123,
  ...
)
```

## Arguments

- object:

  A fitted object. Methods are currently provided for `JoiNMeFit`.

- ...:

  Arguments passed to a class-specific method.

- newdataLong:

  Longitudinal evaluation data; defaults to the fitted data.

- newdataEvent:

  Event evaluation data; defaults to the fitted data.

- time_start:

  Numeric landmark time or vector of landmark times. The default is
  zero.

- time_horizon:

  Numeric horizon time. It may be scalar or have the same length as
  `time_start`. When both `time_horizon` and `Dt` are omitted, the
  fitted maximum time is used.

- Dt:

  Positive horizon width used when `time_horizon` is omitted.

- cause:

  Integer competing-risk cause, with one denoting the primary event
  type.

- n_samples:

  Number of posterior draws used for dynamic prediction.

- seed:

  Integer seed used when posterior draws are subsampled.

## Value

For a `JoiNMeFit`, a data frame with landmark and horizon times,
cumulative/dynamic AUC, and the numbers of cases, controls, and
comparable case-control pairs.

## Details

For case risks \\r_i\\ and control risks \\r_j\\, the estimator is the
proportion of comparable pairs satisfying \\r_i \> r_j\\, with half
credit for ties. This is the empirical cumulative/dynamic AUC. It does
not apply inverse-probability-of-censoring weights; early-censored
subjects are omitted explicitly and the returned counts make the
resulting comparison set clear.
