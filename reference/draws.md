# Posterior draws for JoiNMe objects

Returns posterior draws in `posterior`-compatible formats with
user-facing parameter names.

Compared with
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md),
`draws()` is the convenience layer for downstream posterior workflows.
It is designed for tasks such as `posterior` summarisation, `bayesplot`
visualisation, regex-based variable selection, and any workflow that
expects a standard draws object.

For fitted `JoiNMe` models, known Stan variables are relabelled from raw
fit object to friendly parameter names. For dynamic prediction objects,
stored draw blocks are flattened into a single draws object with
explicit id, scale, marker, and time labels.

Fitted-object draw arrays are cached inside the underlying R6 container
after the first request so later summaries, diagnostics, and plotting
methods can reuse the same renamed draw payload without re-reading the
backend fit.

Use `draws()` when you want:

- one posterior object containing renamed variables,

- [`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)
  and regex-style variable filtering,

- `bayesplot` directly, similar to
  [`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html),

- a standard draws array/matrix/data frame rather than a
  component-specific extraction payload.

- `as.array` is a shorthand for `draws(format = "draws_array")`.

Use
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
instead when you want:

- one model component at a time (`"fixef"`, `"assoc"`, `"gamma_w"`,
  `"basehaz"`, etc.),

- the explicit `term_map` telling you how user-facing labels map back to
  raw Stan variables,

- special structured payloads such as `what = "association_plot"`,

- prediction draw blocks separated by semantic role before flattening.

## Usage

``` r
draws(object, ...)

# S3 method for class 'JoiNMeFit'
draws(
  object,
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  what = c("all", "basehaz", "baseline_hazard"),
  format = c("draws_array", "draws_matrix", "draws_df"),
  ...
)

# S3 method for class 'JoiNMeDynPred'
draws(
  object,
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  format = c("draws_array", "draws_matrix", "draws_df"),
  ...
)

# S3 method for class 'JoiNMeFit'
as.array(x, ...)

# S3 method for class 'JoiNMeDynPred'
as.array(x, ...)
```

## Arguments

- object:

  A supported JoiNMe object.

- ...:

  Unused.

- variables:

  Optional character vector selecting variables after relabelling.

- regex:

  Logical; if `TRUE`, interpret `variables` as regular expressions.

- draws:

  Optional number of posterior draws to retain.

- seed:

  Integer seed used when subsetting draws.

- what:

  For `JoiNMeFit` objects, either `"all"` (default renamed posterior
  variables), `"basehaz"`, or `"baseline_hazard"`.

- format:

  Output format. Supported values are `"draws_array"`, `"draws_matrix"`,
  and `"draws_df"`.

## Value

A posterior draw object in the requested format. The result is a
`posterior`-compatible object with renamed variables, suitable for
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html),
[`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html),
and `bayesplot`-style visualisation.
