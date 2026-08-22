# Extract individual marker weights

Returns the effective marker weights used by marker-aggregated
association terms. The effective value is the declared marker-specific
offset plus the fitted common weight location and marker-specific
departure. When `marker_weights$shared = TRUE`, one shared set is
returned once. Otherwise, each active weighted association term receives
its own set.

Individual weights are intentionally kept out of `summary.JoiNMeFit()`.
The model summary reports only the common location and within-set
spread. This function and [`coef()`](https://rdrr.io/r/stats/coef.html)
return effective weights,
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) returns
their fitted common locations, and
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) returns
only marker-specific departures after subtracting the common location
and declared offset.

## Usage

``` r
marker_weights(object, ...)

# S3 method for class 'JoiNMeFit'
marker_weights(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...)
```

## Arguments

- object:

  A fitted JoiNMe model.

- ...:

  Unused.

- draws:

  Optional number of posterior draws.

- seed:

  Random seed used when subsetting draws.

- digits:

  Number of decimal places used when `summary = TRUE`.

- summary:

  Logical; return posterior summaries when `TRUE`, otherwise a
  draws-by-weight matrix.

## Value

A data frame with one row per marker and weight set when
`summary = TRUE`; otherwise a posterior draw matrix.
