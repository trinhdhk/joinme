# Plot diagnostics, fitted trajectories, and association curves for JoiNMe fits

Plot diagnostics, fitted trajectories, and association curves for JoiNMe
fits

## Usage

``` r
# S3 method for class 'JoiNMeFit'
plot(
  x,
  type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean",
    "running_quantile", "mcmc"),
  pars = NULL,
  regex_pars = NULL,
  draws = 400,
  seed = 1,
  max_vars = 4,
  mcmc_type = c("intervals", "areas", "dens", "dens_overlay", "hist", "trace", "violin",
    "acf", "rhat", "neff"),
  quantile_probs = c(0.1, 0.5, 0.9),
  subject = NULL,
  marker = NA,
  scale = NULL,
  longitudinal_style = c("curves", "heatmap"),
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  threshold = 0.05,
  smooth_trajectory = TRUE,
  smooth_method = c("loess", "spline"),
  smooth_span = 0.3,
  ci_levels = c(0.5, 0.95),
  ci_type = c("ribbon", "line", "both"),
  observed_first = TRUE,
  facet_by = c("marker", "none"),
  facet_scales = "free_y",
  combined = TRUE,
  show_data = TRUE,
  show_observed_line = TRUE,
  observed_style = list(color = "black", shape = 21, size = 2, alpha = 0.6),
  prediction_style = list(color = "steelblue", fill = "steelblue", linewidth = 0.8, alpha
    = 0.2),
  theme_fn = ggplot2::theme_bw,
  palette_marker = NULL,
  association_options = list(),
  ...
)
```

## Arguments

- x:

  A fitted object of class `JoiNMeFit`.

- type:

  Plot type. Diagnostic types are `"rhat"`, `"ess_bulk"`, `"ess_tail"`,
  `"mcse_mean"`, `"mcse_sd"`, `"running_mean"`, and
  `"running_quantile"`. Fitted-data types are `"longitudinal"`,
  `"survival"`, `"cumhaz"`, and `"association"`. `"mcmc"` delegates to
  [`mcmc_plot()`](https://rdrr.io/pkg/joinme/man/mcmc_plot.html) for
  bayesplot-backed posterior displays. The compatibility alias
  `"longitudinal_heatmap"` is treated as
  `type = "longitudinal", longitudinal_style = "heatmap"`.

- pars:

  Optional character vector of parameter names to include for diagnostic
  plots.

- regex_pars:

  Optional regular expression for parameter selection for diagnostic
  plots.

- draws:

  Optional number of posterior draws to subset for diagnostics. For
  fitted-data plotting types (`"longitudinal"`, `"survival"`,
  `"cumhaz"`), this also controls the number of fitted posterior draws
  used to build plot summaries directly from the fitted object.

- seed:

  Random seed for draw subsetting and fitted plotting.

- max_vars:

  Maximum number of parameters for running diagnostics.

- mcmc_type:

  Bayesplot geometry used when `type = "mcmc"`.

- quantile_probs:

  Numeric vector of quantile probabilities for running quantile plots.

- subject:

  Optional vector of subject ids for fitted-data plots.

- marker:

  Optional marker subset for longitudinal and association plots. Use a
  specific marker level to plot only that marker, or `NA` to plot all
  markers.

- scale:

  Longitudinal scale for fitted-data plots. One of `"epred"`,
  `"linpred"`, or `"predict"`.

- longitudinal_style:

  Display style for fitted longitudinal plots. `"curves"` shows the
  existing separated trajectory curves with credible bands. `"heatmap"`
  shows marker-level mean change over time using draw-level aggregation
  across the selected subjects. The heatmap is therefore an alternative
  display for the longitudinal process rather than a separate modelling
  target.

- longitudinal_times:

  Optional numeric vector of times at which fitted longitudinal
  trajectories are evaluated. Times are expressed on the original
  study-time scale. When `NULL`, an evenly spaced sequence spanning the
  selected subjects' observed follow-up is used.

- longitudinal_points:

  Number of evenly spaced evaluation times used when
  `longitudinal_times = NULL`. The default is 80.

- threshold:

  Posterior sign-certainty threshold for
  `longitudinal_style = "heatmap"`. A tile is shown as significant when
  the posterior probability of either a positive or a negative change,
  relative to the earliest plotted time for that marker, is at least
  `1 - threshold`. Tiles that do not meet that criterion remain visible
  at reduced opacity.

- smooth_trajectory:

  Logical; add a secondary smooth of the posterior median trajectory.

- smooth_method:

  Smoothing method, either `"loess"` or `"spline"`.

- smooth_span:

  Numeric span used by the loess smoother.

- ci_levels:

  Numeric credible interval levels.

- ci_type:

  Credible interval display: `"ribbon"`, `"line"`, or `"both"`.

- observed_first:

  Logical; shade the observed-history region for dynamic prediction
  plots.

- facet_by:

  Faceting choice, either `"marker"` or `"none"`.

- facet_scales:

  Scale rule passed to marker facets.

- combined:

  Logical; combine requested panels when possible.

- show_data:

  Logical; display measured longitudinal responses.

- show_observed_line:

  Logical; connect measured responses within marker.

- observed_style:

  Named list controlling measured-response appearance.

- prediction_style:

  Named list controlling posterior trajectory appearance.

- theme_fn:

  Function returning a ggplot2 theme.

- palette_marker:

  Optional marker colour palette.

- association_options:

  Named list of options for `type = "association"`. Supported entries
  are `association_term`, `association_grid`, `association_range`,
  `association_points`, and `association_metric`. When
  `association_grid` is not supplied, `plot.JoiNMeFit()` uses cached
  model-implied raw support from the fitted association channels and
  only falls back to knot support or observed-data heuristics when that
  cache is unavailable. `association_range`, when supplied, overrides
  those default range heuristics and constructs an evenly spaced
  raw-scale grid over the requested interval.
  `association_metric = "hazard"` plots the posterior contribution used
  by the fitted model. For covariance-style channels (`corr`, `vcov`)
  this contribution is zero-referenced at raw value `0` before
  multiplying by \\\alpha\\; `association_metric = "transform"` plots
  the transform \$f(x)\$ alone and therefore omits the
  association-coefficient sign.

- ...:

  Unused.

## Value

A `ggplot` object, a combined plot, or a named list of plots.

## Details

Fitted-data plotting from `JoiNMeFit` uses fitted posterior samples and
does not accept conditioning or future-prediction arguments
(`condition`, `conditioning`, `time_start`, `times`, `time_horizon`, or
`pred_control`). For population-level covariate profiles, call
[`conditional_effects.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/conditional_effects.JoiNMeFit.md).
For subject-specific forecasts conditional on observed marker history,
call [`predict()`](https://rdrr.io/r/stats/predict.html) and then plot
the resulting `JoiNMeDynPred` object. At fitted trajectory times,
covariates other than time are held at their first observed value within
each selected subject-marker combination.
