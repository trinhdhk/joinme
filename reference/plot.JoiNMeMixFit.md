# Plot class-specific latent-progress summaries

Adds four mixture displays to the inherited JoiNMe plotting interface:

- `type = "longitudinal"` with `estimand = "mean_per_class"` evaluates
  the fixed-effects trajectory plus the centre of every class-specific
  random-effect block;

- `estimand = "marginal_per_class"` averages the inverse-link trajectory
  over the estimated within-class random-effect distribution;

- `type = "class_membership"` shows posterior allocation probabilities;

- `type = "covariance_class"` shows class-specific covariance-regression
  curves, with all \\G\\ classes in each panel.

All other plot types are delegated to
[`plot.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/plot.JoinMeFit.md).

## Usage

``` r
# S3 method for class 'JoiNMeMixFit'
plot(
  x,
  type = "longitudinal",
  estimand = NULL,
  draws = 400,
  seed = 1,
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  marginal_samples = 32L,
  ci_levels = c(0.5, 0.95),
  marker = NULL,
  theme_fn = ggplot2::theme_bw,
  ...
)
```

## Arguments

- x:

  A `JoiNMeMixFit` object.

- type:

  Plot type.

- estimand:

  Optional class-specific longitudinal estimand.

- draws:

  Number of posterior draws.

- seed:

  Posterior subsetting seed.

- longitudinal_times:

  Optional trajectory time grid.

- longitudinal_points:

  Number of default trajectory time points.

- marginal_samples:

  Number of within-component Monte Carlo values per posterior draw for
  `estimand = "marginal_per_class"`.

- ci_levels:

  Credible interval levels.

- marker:

  Optional marker subset.

- theme_fn:

  ggplot2 theme function.

- ...:

  Arguments forwarded to inherited plotting methods.

## Value

A ggplot, a combined display, or a named list of ggplots.
