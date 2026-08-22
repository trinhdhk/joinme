# Conditional effects for a fitted joint model

`conditional_effects()` evaluates model-implied changes in the
longitudinal and event processes while holding all predictors not named
in `effects` at explicit conditioning values. Its interface follows
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
where the joint-model structure permits a direct correspondence.

Unlike a univariate regression, a `JoiNMeFit` contains two statistical
processes. The `process` argument therefore selects the longitudinal
expected response, the event relative hazard, or both. Longitudinal
effects may contain only the fixed-effect contribution, add the fitted
marker-level deviation, or average marker-specific predictions over the
selected markers. Subject- and marker-by-subject deviations are excluded
from all three longitudinal estimands unless `reuse_fitted_re = TRUE`
explicitly selects a subject represented in the fitted model. Event
effects describe the proportional-hazards multiplier `exp(W gamma)` (or
its logarithm), with the baseline hazard and longitudinal association
contribution held fixed. Consequently, event contrasts isolate the part
of the event process attributable to covariates in `formulaEvent`.

The returned data follow the `brms` conditional-effects convention. Each
effect table contains `estimate__`, `se__`, `lower__`, `upper__`, and
`cond__`, plus `effect1__` and, for two-predictor effects, `effect2__`.
JoiNMe adds `process__`, `marker__`, `longitudinal_estimand__`, and
`event_type__` where applicable. This permits JoiNMe to delegate the
graphical construction to the tested `brms` plotting method rather than
maintaining a second plotting grammar.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
conditional_effects(
  x,
  effects = NULL,
  conditions = NULL,
  int_conditions = NULL,
  process = c("longitudinal", "event"),
  prob = 0.95,
  robust = TRUE,
  method = c("posterior_epred", "posterior_linpred"),
  longitudinal_estimand = c("population", "marker", "marginal_marker"),
  event_scale = c("hazard_ratio", "log_hazard_ratio"),
  resolution = 100L,
  surface = FALSE,
  markers = NULL,
  draws = NULL,
  summary = TRUE,
  reuse_fitted_re = FALSE,
  seed = 1,
  plot = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted object of class `JoiNMeFit`.

- effects:

  Effects to evaluate. Supply a character vector such as
  `c("time", "time:treatment")` to use the same requests wherever they
  are valid, or a named list with `longitudinal` and `event` components
  to choose process-specific effects. Interactions may contain at most
  two predictors. If `NULL`, all model predictors and fitted two-way
  interactions are used.

- conditions:

  Optional data frame containing values of predictors on which to
  condition. One set of effects is evaluated for every row.
  `conditions$cond__`, when supplied, provides the facet label;
  otherwise row names are used. Tables returned by
  [`make_conditions()`](https://trinhdhk.github.io/joinme/reference/make_conditions.md)
  can be supplied directly.

- int_conditions:

  Optional named list controlling the evaluation values of predictors
  named in `effects`. Each element may be a vector or a function applied
  to the observed predictor. By default, the first numeric predictor
  spans its observed range, factors use all levels, and a second numeric
  predictor is evaluated at its mean and mean plus or minus one standard
  deviation. This mirrors the principal `brms` convention.

- process:

  Character vector selecting `"longitudinal"`, `"event"`, or both. The
  shorthand `"both"` and the default `c("longitudinal", "event")` both
  return the two processes.

- prob:

  Probability covered by the equal-tailed posterior uncertainty
  interval. The default is `0.95`.

- robust:

  Logical. If `TRUE`, posterior medians define `estimate__`; if `FALSE`,
  posterior means are used.

- method:

  Longitudinal posterior scale. Supported values are `"posterior_epred"`
  and `"posterior_linpred"`. The event process is controlled separately
  by `event_scale`.

- longitudinal_estimand:

  Longitudinal quantity to evaluate. `"population"` uses only the
  fixed-effect contribution, corresponding to exclusion of group-level
  effects in `brms`, and averages the selected marker trajectories
  within each posterior draw. `"marker"` adds the fitted marker-level
  deviation and retains a separate trajectory for every selected marker.
  `"marginal_marker"` also adds the marker-level deviation but then
  takes an equally weighted mean across selected markers within every
  draw. Marker-specific inverse links are applied before either marker
  average. Subject and marker-by-subject deviations are excluded unless
  `reuse_fitted_re = TRUE`.

- event_scale:

  Event-process estimand. `"hazard_ratio"` returns the
  covariate-specific proportional-hazards multiplier `exp(W gamma)`;
  `"log_hazard_ratio"` returns `W gamma`.

- resolution:

  Number of support points for the first continuous predictor. For a
  two-dimensional surface it is used on both axes.

- surface:

  Logical. If `TRUE`, two continuous predictors are evaluated over a
  rectangular surface. If `FALSE`, the second continuous predictor is
  represented by a small set of conditioning curves.

- markers:

  Optional character vector restricting longitudinal results to selected
  marker levels. By default all fitted markers are used. For
  `longitudinal_estimand = "population"` or `"marginal_marker"`, this
  argument defines the markers entering the equally weighted posterior
  average.

- draws:

  Optional positive integer limiting the posterior draws used in the
  calculation. `NULL` uses all available draws.

- summary:

  Logical. If `TRUE`, return posterior centres, standard errors, and
  credible intervals. If `FALSE`, return a tidy draw-level data frame.
  Each row contains one posterior value, its `.draw` index, its
  `.value`, and the complete condition and estimand description.
  Draw-level results are returned without plotting.

- reuse_fitted_re:

  Logical. If `TRUE`, `conditions` must contain the fitted subject
  identifier. Subject and marker-by-subject posterior effects for those
  identifiers are added draw by draw, avoiding a separate dynamic
  random-effect fit. All identifiers must have occurred during fitting.

- seed:

  Integer seed used when posterior draws are subsampled.

- plot:

  Logical. If `TRUE` and `summary = TRUE`, construct and display
  conditional-effects plots. If `FALSE`, return the underlying
  `JoiNMeConditionalEffects` object. `summary = FALSE` always returns
  the draw-level object without plotting.

- ...:

  Additional arguments passed to the plotting method when `plot = TRUE`,
  for example `points`, `rug`, `stype`, `theme`, or the conditional-plot
  arrangement arguments; see
  [`plot.JoiNMeConditionalEffects()`](https://trinhdhk.github.io/joinme/reference/plot.JoiNMeConditionalEffects.md)
  and
  [`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html).

## Value

If `summary = TRUE` and `plot = FALSE`, a `JoiNMeConditionalEffects`
object: a named list with one component per requested process. Each
process component is a `brms_conditional_effects`-compatible named list
containing one data frame per effect. Longitudinal tables identify their
estimand in `longitudinal_estimand__`. With `summary = FALSE`, each
effect is a tidy data frame with one row per posterior draw and
estimand, identified by `.draw`, `.value`, and `estimand__`. If
`summary = TRUE` and `plot = TRUE`, the plotting method invisibly
returns either the similarly nested named list of `ggplot` objects or an
arranged `patchwork` display.

## Details

Numeric predictors not named in an effect are held at their observed
mean; factors are held at their first level. Values supplied in
`conditions` override those defaults. Group-level subject deviations are
excluded, as in the default `re_formula = NA` behavior of
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html).
Marker-level deviations are excluded for the `"population"` estimand,
retained separately for the `"marker"` estimand, and integrated by
finite averaging for the `"marginal_marker"` estimand. Population and
marker-marginal results each contain one trajectory per condition
because their selected marker trajectories are averaged draw by draw.
When fitted random effects are reused, the selected subject and
marker-by-subject deviations enter all three estimands.

The event result is a relative-hazard effect rather than a dynamic
survival prediction. Dynamic survival probabilities depend on a
subject's observed marker history and conditioning time and therefore
remain the estimand of
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md).
Keeping these two estimands separate prevents a table of baseline
covariate effects from being mislabeled as subject-specific dynamic
prediction.

## See also

[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html),
[`make_conditions()`](https://trinhdhk.github.io/joinme/reference/make_conditions.md),
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Default longitudinal and event effects, returned as plot-ready data.
cond_eff <- conditional_effects(fit, plot = FALSE)

# Compare treatment profiles for a longitudinal time effect and an event
# age effect, using labels generated by brms::make_conditions().
profiles <- make_conditions(fit$dataEvent, vars = "treatment")
cond_eff <- conditional_effects(
  fit,
  effects = list(longitudinal = "time", event = "age"),
  conditions = profiles,
  process = c("longitudinal", "event"),
  plot = FALSE
)
plot(cond_eff, arrange = "grid", ncol = 2, guides = "collect")

# Request marker-specific longitudinal expected responses.
conditional_effects(
  fit,
  effects = "time:treatment",
  process = "longitudinal",
  longitudinal_estimand = "marker",
  int_conditions = list(treatment = c("control", "active"))
)

# Average marker-specific expected responses draw by draw.
conditional_effects(
  fit,
  effects = "time",
  process = "longitudinal",
  longitudinal_estimand = "marginal_marker",
  markers = c("marker_1", "marker_2"),
  plot = FALSE
)
} # }
```
