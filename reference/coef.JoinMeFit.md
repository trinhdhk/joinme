# Combined posterior coefficients for fitted JoiNMe models

Returns posterior coefficients on the scale used by each model
component.

The guiding rule is simple: for every coefficient carried by a
group-specific model matrix, the returned value is the sum of the
population-level contribution and the matching group-level deviation.
When no group-level deviation exists, the returned coefficient is the
population-level coefficient itself.

This mirrors the interpretation used in multilevel modelling: a
subject-specific or marker-specific coefficient is the coefficient that
would multiply the corresponding column of the model matrix for that
unit.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
coef(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- draws:

  Optional number of posterior draws to retain.

- seed:

  Integer seed used when subsetting posterior draws.

- digits:

  Number of digits used when `summary = TRUE`.

- summary:

  Logical. If `TRUE`, return posterior summaries. If `FALSE`, return
  posterior draw-level extractions.

- ...:

  Unused.

## Value

When `summary = FALSE`, a nested list of draw-level data frames for the
longitudinal, event, distributional, and covariance-regression parts of
the model. When `summary = TRUE`, the same structure is returned after
summarising each coefficient with posterior means, posterior
uncertainty, interval estimates, and MCMC diagnostics.
