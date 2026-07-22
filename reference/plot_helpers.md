# Convenience plot helpers for JoiNMe objects

Convenience wrappers around the main
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
JoiNMe fit and dynamic prediction objects.

[`longitudinal_plot()`](https://trinhdhk.github.io/joinme/reference/longitudinal_plot.html),
[`survival_plot()`](https://trinhdhk.github.io/joinme/reference/plot_helpers.html),
and
[`cumhaz_plot()`](https://trinhdhk.github.io/joinme/reference/cumhaz_plot.html)
route to the corresponding fitted or predicted trajectory plot.

[`association_plot()`](https://trinhdhk.github.io/joinme/reference/association_plot.html)
exposes the fitted association-curve display for `JoiNMeFit` objects.

[`diagnostic_plot()`](https://trinhdhk.github.io/joinme/reference/diagnostic_plot.html)
exposes the scalar and running sampler diagnostics for `JoiNMeFit`
objects.

[`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html)
applies bayesplot MCMC geometries to the renamed posterior draws
returned by
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md), which
means the displayed parameter labels are the same user-facing names used
by summaries and diagnostics.

## Usage

``` r
longitudinal_plot(object, ...)

survival_plot(object, ...)

cumhaz_plot(object, ...)

association_plot(object, ...)

diagnostic_plot(
  object,
  type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean",
    "running_quantile"),
  ...
)

mcmc_plot(
  object,
  pars = NA,
  type = c("intervals", "areas", "dens", "dens_overlay", "hist", "trace", "violin",
    "acf", "rhat", "neff"),
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  draws = NULL,
  seed = 1,
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit`, `JoiNMeFit`, `JoiNMeDynPred`, or `PredJoiNMeFit` object
  as appropriate for the requested helper.

- type:

  The requested diagnostic or bayesplot geometry.

- pars:

  Deprecated alias of `variable` for
  [`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html).

- variable:

  Optional renamed posterior variable names for
  [`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html).

- regex:

  Logical; whether `variable` should be treated as a regular expression.

- fixed:

  Deprecated logical alias controlling exact versus regular-expression
  matching when `pars` is supplied.

- draws:

  Optional number of posterior draws to keep.

- seed:

  Integer seed used when subsetting posterior draws.

- ...:

  Additional arguments forwarded to the corresponding
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
  `bayesplot::mcmc_*()` implementation.

## Value

For the plot helpers, a `ggplot` object, a combined plot object, or a
named list of plots depending on the requested object and plot family.
