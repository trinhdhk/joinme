# Plot posterior MCMC summaries using bayesplot

Provides a JoiNMe-friendly wrapper around `bayesplot::mcmc_*` functions.
The posterior draws are first relabelled with user-facing parameter
names via
[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md),
after which the selected bayesplot geometry is applied.

## Usage

``` r
mcmc_plot(
  object,
  pars = NA,
  type = c("intervals", "areas", "dens", "dens_overlay", "hist", "trace", "violin",
    "acf", "rhat", "neff"),
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  draws = NULL,
  seed = .Random.seed[[1]],
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit` or `JoiNMeDynPred` object.

- pars:

  Deprecated alias of `variable`.

- type:

  Plot type. Supported values are `"intervals"`, `"areas"`, `"dens"`,
  `"dens_overlay"`, `"hist"`, `"trace"`, `"violin"`, `"acf"`, `"rhat"`,
  and `"neff"`.

- variable:

  Optional variable names after relabelling.

- regex:

  Logical; treat `variable` as a regular expression.

- fixed:

  Deprecated logical alias controlling exact versus regex matching when
  `pars` is supplied.

- draws:

  Optional number of posterior draws to keep.

- seed:

  Integer seed used when subsetting draws.

- ...:

  Additional arguments passed to the selected `bayesplot` function.

## Value

A `ggplot` object.
