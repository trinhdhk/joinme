# Conditional contrasts for fitted JoiNMe models

`conditional_contrast()` compares two named covariate profiles within
every posterior draw and only then summarises the resulting contrast.
This ordering retains posterior dependence between the two
counterfactual predictions and gives the appropriate uncertainty for
their paired difference or hazard ratio.

`groupA` and `groupB` describe the covariate values that distinguish the
two profiles. Other model predictors are taken from each row of
`conditions`, when supplied, or from reference values of the fitting
data. Group values take precedence when a variable also appears in
`conditions`. A named list is useful when values have different types; a
named atomic vector is accepted for concise specifications such as
`c(treatment = "active", sex = "female")`.

For the longitudinal process, the reported contrast is group A minus
group B on the scale selected by `method`. For the event process,
`event_scale = "hazard_ratio"` reports the direct covariate hazard ratio
of group A relative to group B, whereas `"log_hazard_ratio"` reports its
logarithm. As in
[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md),
the event calculation concerns the covariate component of
`formulaEvent`; the baseline hazard and longitudinal association
contribution are held common between profiles.

## Usage

``` r
conditional_contrast(x, ...)

# S3 method for class 'JoiNMeFit'
conditional_contrast(
  x,
  groupA,
  groupB,
  conditions = NULL,
  process = c("longitudinal", "event"),
  prob = 0.95,
  robust = TRUE,
  method = c("posterior_epred", "posterior_linpred"),
  longitudinal_estimand = c("population", "marker", "marginal_marker"),
  event_scale = c("hazard_ratio", "log_hazard_ratio"),
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

  A fitted object inheriting from `JoiNMeFit`.

- ...:

  Additional arguments passed to
  [`plot.JoiNMeConditionalContrasts()`](https://trinhdhk.github.io/joinme/reference/plot.JoiNMeConditionalContrasts.md)
  when `plot = TRUE`.

- groupA:

  Named atomic vector or named list containing one scalar value for
  every covariate fixed to define group A.

- groupB:

  Named atomic vector or named list containing one scalar value for
  every covariate fixed to define group B.

- conditions:

  Optional data frame containing the common predictor values at which
  the contrast is evaluated. There is one posterior contrast for every
  row. A `cond__` column supplies display labels; otherwise informative
  row names or numbered labels are used. Values not supplied by either a
  group or a condition row are derived from the fitting data in the same
  way as
  [`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md).

- process:

  Character vector selecting `"longitudinal"`, `"event"`, or both.
  `"survival"` is accepted as an alias for `"event"`, and `"both"`
  selects both fitted processes.

- prob:

  Probability covered by the equal-tailed posterior uncertainty
  interval. The default is `0.95`.

- robust:

  Logical. If `TRUE`, posterior medians define `estimate__`; if `FALSE`,
  posterior means are used.

- method:

  Longitudinal posterior scale. `"posterior_epred"` compares expected
  responses after applying the marker-specific inverse link;
  `"posterior_linpred"` compares longitudinal linear predictors.

- longitudinal_estimand:

  Longitudinal quantity entering the comparison. `"population"` uses
  population coefficients only and averages the selected marker
  trajectories within each draw. `"marker"` also includes the fitted
  marker-level deviation and retains a contrast for each selected
  marker. `"marginal_marker"` includes that deviation and averages the
  selected marker-specific predictions within each draw. Subject and
  marker-by-subject deviations are excluded unless
  `reuse_fitted_re = TRUE` selects a fitted subject.

- event_scale:

  Event-process contrast. `"hazard_ratio"` reports `exp(eta_A - eta_B)`
  and `"log_hazard_ratio"` reports `eta_A - eta_B`, where `eta` is the
  direct event-regression linear predictor.

- markers:

  Optional character vector restricting longitudinal marker levels. With
  `longitudinal_estimand = "population"` or `"marginal_marker"`, these
  are the marker trajectories averaged within each posterior draw.

- draws:

  Optional positive integer limiting the posterior draws used in the
  paired calculation. `NULL` uses all available draws.

- summary:

  Logical. If `TRUE`, return posterior centres, standard errors, and
  credible intervals. If `FALSE`, each process contains a tidy data
  frame with one row per posterior draw and estimand. The `.draw` and
  `.value` columns identify the retained draw and its value. Draw-level
  results are returned without plotting.

- reuse_fitted_re:

  Logical. If `TRUE`, `conditions` must contain a fitted subject
  identifier. The paired comparison then includes that subject's
  posterior subject and marker-by-subject effects without fitting new
  random effects. The identifier cannot be changed between `groupA` and
  `groupB`.

- seed:

  Integer seed used when posterior draws are subsampled.

- plot:

  Logical. If `TRUE` and `summary = TRUE`, construct and display
  contrast plots. If `FALSE`, return the `JoiNMeConditionalContrasts`
  object. `summary = FALSE` always returns draw-level results without
  plotting.

## Value

With `summary = TRUE` and `plot = FALSE`, a `JoiNMeConditionalContrasts`
object containing one data frame per requested process. Every table
includes `estimate__`, `se__`, `lower__`, and `upper__`, together with
condition, group, marker, event-type, and estimand labels where
applicable. With `summary = FALSE`, every process is instead a tidy
draw-level data frame containing `.draw`, `.value`, `estimand__`, and
the applicable scientific descriptors. With `summary = TRUE` and
`plot = TRUE`, the plotting method invisibly returns either the named
list of `ggplot` objects or an arranged `patchwork` display.

## Details

Factor levels, transformations, spline bases, time scaling, and
coefficient ordering are recovered from the fitted model. The two
profile designs are evaluated in one posterior calculation so their
columns refer to precisely the same MCMC draws. The contrast is then
formed draw by draw. Credible intervals therefore describe the posterior
distribution of the contrast, rather than a subtraction of separately
summarised intervals.

If `conditions` is omitted, numeric predictors use their observed mean
and categorical predictors use their first fitted level. Multiple
condition rows may describe distinct profiles or a trajectory over a
continuous variable such as follow-up time. The plotting method can
select such a variable with `condition_variable`.

## See also

[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md),
[`make_conditions()`](https://trinhdhk.github.io/joinme/reference/make_conditions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
treatment_contrast <- conditional_contrast(
  fit,
  groupA = c(treatment = "active"),
  groupB = c(treatment = "control"),
  process = c("longitudinal", "event"),
  plot = FALSE
)
print(treatment_contrast)
plot(treatment_contrast, arrange = "row", guides = "collect")

time_profiles <- data.frame(
  time = seq(0, 24, length.out = 50),
  cond__ = paste0("month ", seq(0, 24, length.out = 50))
)
treatment_over_time <- conditional_contrast(
  fit,
  groupA = list(treatment = "active", sex = "female"),
  groupB = list(treatment = "control", sex = "female"),
  conditions = time_profiles,
  process = "longitudinal",
  method = "posterior_epred",
  longitudinal_estimand = "marginal_marker",
  plot = FALSE
)
plot(treatment_over_time, condition_variable = "time")
} # }
```
