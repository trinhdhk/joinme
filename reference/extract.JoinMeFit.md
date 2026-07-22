# Extract posterior draws from a fitted JoiNMe model

Extracts one fitted-model component at a time and returns both the
draw-level values and the mapping that produced them.

This differs from
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md) in two
important ways:

- [`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
  keeps the request scoped to one semantic component such as fixed
  effects, survival coefficients, association terms, distributional
  terms, or likelihood-scale parameters.

- [`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
  returns a structured list with `draws` plus `term_map` (and for
  `what = "association_plot"`, additional plotting support data) instead
  of a single `posterior` draws object.

In short, use
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
when you need component-aware extraction and use
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md) when
you need one renamed posterior object for general downstream analysis.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
extract(
  object,
  what = c("fixef", "gamma_w", "basehaz", "baseline_hazard", "assoc", "association_plot",
    "distributional", "distributional_regression", "likelihood_scale", "raw"),
  term = NULL,
  variable = NULL,
  draws = NULL,
  seed = 1,
  keep_chains = TRUE,
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- what:

  Character component selector. One of `"fixef"`, `"gamma_w"`,
  `"basehaz"`, `"baseline_hazard"`, `"assoc"`, `"association_plot"`,
  `"distributional"`, `"distributional_regression"`,
  `"likelihood_scale"`, or `"raw"`.

- term:

  Optional character vector of friendly term names (summary-style) to
  subset extracted columns.

- variable:

  Optional character vector of raw Stan variable names. This is used
  directly when `what = "raw"` and can also further filter mapped
  outputs.

- draws:

  Optional number of posterior draws to keep per chain.

- seed:

  Integer seed used when subsetting draws.

- keep_chains:

  Logical; if TRUE, return draws with chains in a separate dimension
  (iteration x chain x term). If FALSE, return a flattened draws-by-term
  matrix.

## Value

A list with fields:

- `draws`: numeric array when `keep_chains = TRUE` (iteration x chain x
  term), otherwise a numeric matrix (rows = draws, cols = requested
  terms). For `what = "association_plot"`, this is a named list of
  compact draw matrices keyed by association term.

- `term_map`: data.frame mapping `term` to Stan `variable`

- `support`: for `what = "association_plot"`, cached model-implied raw
  support ranges used by association plotting.

## See also

[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md) for a
higher-level interface that returns a single `posterior`
