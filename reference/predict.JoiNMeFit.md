# Dynamic Prediction for Joint Models

Performs dynamic prediction for a new subject using a fitted `JoiNMeFit`
joint model. It calculates the posterior predictive distribution of
longitudinal trajectories and survival probabilities conditional on the
subject's observed history up to a specific time point.

Convenience wrapper around the
[predict](https://rdrr.io/r/stats/predict.html) method returning the
posterior linear predictor (linpred scale).

Convenience wrapper around the
[predict](https://rdrr.io/r/stats/predict.html) method returning the
posterior expected predictor (epred scale).

Convenience wrapper around the
[predict](https://rdrr.io/r/stats/predict.html) method returning
posterior predictive draws (includes observation noise).

## Usage

``` r
# S3 method for class 'JoiNMeFit'
predict(
  object,
  newdataLong,
  newdataEvent = NULL,
  process = c("longitudinal", "event"),
  pred_type = c("per_marker_id", "marginal_marker", "marginal_id", "marker_subject",
    "subject_marker"),
  scale = c("epred", "linpred", "predict"),
  times = NULL,
  time_start = NULL,
  time_horizon = NULL,
  tmax = NULL,
  ci_levels = c(0.5, 0.95),
  control = list(),
  reuse_fitted_re = FALSE,
  seed = .Random.seed[[1]],
  ...
)

# S3 method for class 'JoiNMeFit'
posterior_linpred(object, reuse_fitted_re = FALSE, ...)

# S3 method for class 'JoiNMeFit'
posterior_epred(object, reuse_fitted_re = FALSE, ...)

# S3 method for class 'JoiNMeFit'
posterior_predict(object, reuse_fitted_re = FALSE, ...)
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

- pred_type:

  Character. Type of longitudinal predictions:

  - "per_marker_id" (default): Subject-specific predictions for each
    marker.

  - "marginal_marker": Average across subjects for each marker
    (population-level by marker).

  - "marginal_id": Average across markers for each subject
    (subject-specific average).

  - "marker_subject": Alias for "marginal_marker".

  - "subject_marker": Alias for "marginal_id".

- scale:

  Character vector. Longitudinal prediction scale(s):

  - "epred": expected response (inverse-link),

  - "linpred": linear predictor,

  - "predict": predictive draw (includes noise). Defaults to all scales
    when not explicitly set. For `epred`/`predict`, family-specific
    inverse-links configured in `families` (via `jm_family(...)`) are
    respected.

- times:

  Numeric vector of times at which to predict the longitudinal and
  survival trajectories. Can also be a named list of numeric vectors
  (one per subject id). If NULL, a grid from `time_start` to
  `time_start + time_horizon` is generated (and truncated to training
  support when needed). If fewer than 50 points are supplied for a
  subject, a warning is emitted when `control$n_times = 50`.

- time_start:

  Numeric scalar or character. The conditioning time (last observation
  time). If numeric, a single value is reused for all subjects, or a
  named vector supplies subject-specific values. If character, it is
  interpreted as a column name in `newdataEvent` holding
  subject-specific conditioning times. If NULL, defaults to the maximum
  observed time in `newdataLong` for each subject; for delayed-entry
  counting-process data with strictly positive entry times, it defaults
  to the subject entry time.

- time_horizon:

  Numeric scalar. Prediction horizon (in time units) used when `times`
  is NULL. The default is `tmax`.

- tmax:

  Numeric scalar. The maximum time used for scaling during model
  fitting. If this is not provided and cannot be inferred from the
  object, predictions will be on the wrong time scale.

- ci_levels:

  Numeric vector of credible interval levels for plotting. Must be
  strictly between 0 and 1.

- control:

  Named list for prediction configuration.

  - cmdstanr::model\$sample() arguments (e.g., `chains`,
    `parallel_chains`, `iter_warmup`, `iter_sampling`, `seed`,
    `refresh`, `adapt_delta`). Values override defaults except `data`.

  - engine: "cmdstanr" or "rstan". Defaults to
    options(stan_preferred_engine).

  - n_samples: integer; number of posterior parameter draws extracted
    from the fitted object before dynamic prediction. Default 200.

  - n_times: integer; number of points in the prediction time grid when
    `times` is `NULL`. Default 50.

  - n_pred_draws: integer; number of prediction draws produced by the
    dynpred object. Default to n_samples. Changing this can destabilise
    the results. This is independent of `n_samples` (posterior parameter
    draw extraction count).

  - threads_per_chain: integer; the threaded dynpred Stan program is
    always used. `threads_per_chain = 1` keeps execution serial while
    preserving the thread-capable kernel. For engine = "rstan",
    threading uses options(stan.thread = threads_per_chain).

  - grainsize: integer; reduce_sum grainsize for threaded prediction.
    Defaults to the full prediction draw count when
    `threads_per_chain = 1`, and to
    `max(1, ceiling(n_pred_draws/(4*threads_per_chain*chains)))`
    otherwise.

  - quadrature_nodes: optional positive integer target for total
    quadrature points during dynamic prediction. Allowed values are
    exactly `7`, `15`, `31`, `41`, `51`, and `61`. Only the node count
    is passed to Stan; GK nodes/weights are fixed in the Stan code.

  - progress: logical; show sampling progress bar (default TRUE).

  - initialisation fallback (cmdstanr): if a subject-level dynamic
    prediction chain fails to initialise, prediction automatically
    retries that subject with a narrower random init range
    (`init = 0.1`).

- reuse_fitted_re:

  Logical. If `TRUE`, every identifier in the supplied data must have
  occurred during fitting. Prediction reuses that subject's paired
  posterior random effects and covariance draws, without estimating a
  second set of random effects. The default, `FALSE`, conditions newly
  sampled effects on the supplied longitudinal history for dynamic
  prediction.

- seed:

  Integer. Random seed for reproducibility of random effect sampling.

- ...:

  Additional arguments passed to
  [predict](https://rdrr.io/r/stats/predict.html).

## Value

A list with two components:

- longitudinal:

  A data.frame containing longitudinal predictions (Mean, Median, SD,
  95% CrI) for each time point in `times`. Columns: id, time, marker,
  Estimate, Median, Est.Error, L95, U95.

- longitudinal_fitted:

  A data.frame containing fitted values for observed history on the
  linpred and epred scales. Includes a `scale` column.

- survival:

  A data.frame containing survival probabilities (S(t\|time_start)) for
  each time point in `times`. Columns: id, time, Survival, Median,
  Est.Error, L95, U95.

- cumhaz:

  A data.frame containing conditional cumulative hazards
  H(t\|time_start). Columns: id, time, Cumhaz, Median, Est.Error, L95,
  U95.

- draws:

  A list containing raw posterior draws (`longitudinal`,
  `longitudinal_fitted`, `survival`, `cumhaz`) and reconstructed
  subject-level random effects (`random_effects_id`,
  `random_effects_marker_id`, including marker-by-id covariance draws
  when id-dependent covariance is active). Predictions from
  [`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md)
  additionally retain conditional allocation draws in `posterior_class`.

A `JoiNMeDynPred` object with `metadata$scale = "linpred"` and
`metadata$scales = "linpred"`.

A `JoiNMeDynPred` object with `metadata$scale = "epred"` and
`metadata$scales = "epred"`.

A `JoiNMeDynPred` object with `metadata$scale = "predict"` and
`metadata$scales = "predict"`.

## Details

`predict.JoiNMeFit()` does not accept a `condition` argument. It
estimates subject-specific future trajectories and survival conditional
on the longitudinal histories supplied in `newdataLong`. Use
[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md)
with
[`make_conditions()`](https://trinhdhk.github.io/joinme/reference/make_conditions.md)
when the target is a population-level comparison across named covariate
profiles.

The function uses a Bayesian approach in two draw layers:

1.  Extract posterior parameter draws from the fitted object
    (`n_samples`).

2.  Re-index those draws to the dynpred draw count
    (`control$n_pred_draws`), allowing independent control of prediction
    Monte Carlo size.

3.  For each dynpred draw, sample subject-specific random effects
    *conditional* on observed history in `newdataLong`.

4.  Evaluate longitudinal and survival quantities on requested future
    grids.

5.  Pool draws to form marginal predictive summaries and intervals.

Dynamic prediction passes the fitted, draw-specific transform ordinates
to Stan. Ordered piecewise-linear associations therefore use the same
posterior curve as the fitted event model and never reconstruct a curve
from earlier user-supplied `y` values.

The baseline-hazard basis is likewise inherited from fitting. Prediction
reuses the retained B-spline or natural-spline object, or the retained
baseline formula, together with the original column-centring constants.
This guarantees that the quadrature arrays have the same number and
meaning of baseline-hazard columns as the fitted Stan parameters.
