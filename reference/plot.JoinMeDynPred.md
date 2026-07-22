# Enhanced Plot Method for Dynamic Prediction from Joint Models

Produces publication-ready ggplot2 visualisations of dynamic predictions
from a joint model with flexible customisation of longitudinal
trajectories and conditional survival curves.

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
plot(
  x,
  type = c("cumhaz", "longitudinal", "survival"),
  subject = NULL,
  marker = NA,
  trajectory_type = c("id_marker", "marker_pop", "overall_pop", "all_types"),
  scale = NULL,
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
  ...
)
```

## Arguments

- x:

  An object of class `JoiNMeDynPred` returned by the
  [predict](https://rdrr.io/r/stats/predict.html) method.

- type:

  Character vector indicating what to plot: "cumhaz" (default),
  "longitudinal", or "survival". Multiple values are allowed.

- subject:

  Integer/Character vector. Which subject(s) to plot. If NULL, plots all
  subjects in separate panels or returns a list.

- marker:

  Optional marker subset for longitudinal plots. Use a specific marker
  level to plot only that marker, or `NA` to plot all markers.

- trajectory_type:

  Character. Which trajectory to show: "id_marker" (default),
  "marker_pop" (marker-level average), "overall_pop" (overall mean), or
  "all_types" for all three.

- scale:

  Optional longitudinal prediction scale to plot. Must match available
  scales in the prediction object ("epred", "linpred", "predict"). If
  NULL, defaults to "epred" when available, otherwise the first
  available scale.

- smooth_trajectory:

  Logical. Whether to smooth trajectories with splines (default TRUE).
  Only applicable if sufficient observations.

- smooth_method:

  Character. Smoothing method: "loess" or "spline".

- smooth_span:

  Numeric in (0,1). Controls span/flexibility of LOESS smoothing
  (default 0.3).

- ci_levels:

  Numeric vector. Credible interval levels to display (default c(0.5,
  0.95) for 50% and 95%).

- ci_type:

  Character. How to display uncertainty: "ribbon" (shaded area), "line"
  (quantile lines), or "both".

- observed_first:

  Logical. Whether to visually separate observed history (shaded
  background) from prediction region (default TRUE).

- facet_by:

  Character. "marker" to facet by marker type, "none" for single plot.

- facet_scales:

  Character. "fixed" or "free_y" for longitudinal faceting.

- combined:

  Logical. If TRUE and multiple plot types are requested, combine them
  into a single layout (patchwork/cowplot).

  - For a single subject, returns one combined plot.

  - For multiple subjects, returns a named list of combined plots (one
    per subject).

  - If a combiner package is unavailable, falls back to the uncombined
    subject/outcome structure for that subject.

- show_data:

  Logical. Whether to overlay observed data points (default TRUE).

- show_observed_line:

  Logical. Whether to connect observed points with lines (default TRUE).

- observed_style:

  List with elements "color", "shape", "size", "alpha" for appearance of
  observed data points.

- prediction_style:

  List with elements "color", "fill", "linewidth", "alpha" for
  prediction line and ribbons.

- theme_fn:

  ggplot2 theme function (default
  [`ggplot2::theme_minimal`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).

- palette_marker:

  Character vector of colors for different markers, or function like
  [`ggplot2::scale_color_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html).

- ...:

  Additional arguments (unused).

## Value

If a single subject and outcome requested: a `ggplot` object. If
multiple subjects and `combined = TRUE`: a named list where each element
is a per-subject combined plot (patchwork/cowplot) when backend support
exists, otherwise the corresponding uncombined per-subject list. If
multiple subjects and `combined = FALSE`: a named list grouped by
subject, each containing requested outcome plots. If multiple outcomes:
list with "longitudinal", "survival", and/or "cumhaz" elements, or a
single combined plot if `combined = TRUE`.

## Details

**Longitudinal Plots:**

- The longitudinal scale is selected via `.arg scale` (or inferred when
  NULL).

- On \\predict\\ scale: observed values are shown for times in interval
  `[0, time_start]` and predictions are drawn for
  `[time_start, time_horizon]`.

- On \\linpred\\ or \\epred\\ scale: predictions are shown for
  `[0, time_horizon]`, with fitted values providing the pre-time_start
  segment when available.

- Optional shaded region distinguishing observed from prediction
  periods.

- Predicted trajectories with credible bands.

- Multiple CI levels with decreasing alpha for visual hierarchy.

- Optional smoothing for smoother appearance.

- Support for multiple markers with automatic faceting.

**Survival / Cumulative Hazard Plots:**

- Conditional survival probability S(t \| T_cond).

- Conditional cumulative hazard H(t \| T_cond) (default).

- Credible bands with multiple levels.
