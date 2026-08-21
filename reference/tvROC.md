# Time-dependent ROC curves and areas

Re-exports the `JMbayes2` generics
[`JMbayes2::tvROC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html)
and
[`JMbayes2::tvAUC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html).
JoiNMe supplies methods for fitted `JoiNMeFit` objects, while ROC
objects retain class `"tvROC"` and can therefore use the plotting,
printing, and AUC methods supplied by `JMbayes2`.

`tvROC.JoiNMeFit()` estimates a cumulative/dynamic ROC curve using
longitudinal information available through a landmark time. Cases
experience the requested cause in the interval from the landmark to the
horizon; controls remain event-free beyond the horizon. Censoring is
handled by model-based expected status or Kaplan–Meier
inverse-probability weights.

Dynamic risks are calculated by
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md).
They therefore include the fitted baseline hazard, survival covariates,
longitudinal trajectories, marker weights, and every identity,
functional, monotone-spline, or ordered piecewise-linear association.

## Usage

``` r
tvROC(object, newdata, Tstart, ...)

tvAUC(object, newdata, Tstart, ...)

# S3 method for class 'JoiNMeFit'
tvROC(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  time_start = 0,
  time_horizon = NULL,
  Dt = NULL,
  cause = 1,
  n_samples = 200,
  seed = 123,
  type_weights = c("model-based", "IPCW"),
  ...
)

# S3 method for class 'JoiNMeFit'
tvAUC(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  time_start = 0,
  time_horizon = NULL,
  Dt = NULL,
  cause = 1,
  n_samples = 200,
  seed = 123,
  type_weights = c("model-based", "IPCW"),
  ...
)
```

## Arguments

- object:

  A fitted `JoiNMeFit` object.

- ...:

  Additional arguments passed to
  [`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md).

- newdataLong:

  Longitudinal evaluation data; defaults to the fitted data.

- newdataEvent:

  Event-process evaluation data; defaults to the fitted data.

- time_start:

  Numeric landmark time or vector of landmark times. The default is
  zero.

- time_horizon:

  Numeric horizon. It may be scalar or aligned with `time_start`. When
  this and `Dt` are omitted, the fitted maximum follow-up time is used.

- Dt:

  Positive prediction-window width used when `time_horizon` is omitted.

- cause:

  Positive integer identifying the event cause of interest. Other causes
  are treated as censoring for cause-specific discrimination.

- n_samples:

  Positive integer number of posterior draws used for dynamic
  prediction.

- seed:

  Integer seed used when posterior draws are subsampled.

- type_weights:

  Censoring treatment: `"model-based"` imputes the expected case/control
  status of subjects censored before the horizon; `"IPCW"` applies
  inverse Kaplan–Meier censoring weights to observable cases and
  controls.

## Value

`tvROC.JoiNMeFit()` returns a `"tvROC"` object, or a
`"tvROC_JoiNMeFit_list"` for multiple landmarks. `tvAUC.JoiNMeFit()`
returns the standard `"tvAUC"` object for one landmark and a data frame
of areas for multiple landmarks.

## Details

For a threshold \\c\\, a subject is classified as a case when their
predicted conditional survival probability is below \\c\\. Sensitivity
and one minus specificity are evaluated at 101 thresholds from zero to
one, following
[`JMbayes2::tvROC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html).
Posterior-draw curves are retained in `tp` and `fp`; `TP` and `FP` are
calculated from posterior mean risks.

With model-based weighting, an observed case has case weight one, a
subject known to be event-free beyond the horizon has case weight zero,
and a subject censored within the prediction window contributes their
model-based conditional event probability. IPCW uses the Kaplan–Meier
estimator of the censoring survival distribution conditional on
remaining under observation at the landmark.

A scalar landmark returns a standard `"tvROC"` object compatible with
[`JMbayes2::tvAUC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Multiple
landmarks return a `"tvROC_JoiNMeFit_list"` containing one compatible
curve per landmark–horizon pair.

## References

Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
curves. *Biometrics*, 61, 92–105.

Rizopoulos D (2011). Dynamic predictions and prospective accuracy in
joint models. *Biometrics*, 67, 819–829.

## See also

[`JMbayes2::tvROC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html),
[`JMbayes2::tvAUC()`](https://drizopoulos.github.io/JMbayes2/reference/accuracy.html),
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md),
[`concordance.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/concordance.JoiNMeFit.md)
