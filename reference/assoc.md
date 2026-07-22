# Build posterior association effects for a fitted JoiNMe model

Reconstructs the posterior association effects that enter the survival
linear predictor.

For weighted current-value and current-slope channels (`cv_total`,
`cs_total`, `cv_marker`, `cs_marker`), the returned effect is the
draw-wise product of the association coefficient and the marker weight,
divided by the number of markers to match the scale used in the fitted
hazard contribution.

For scalar channels (`cv_mean`, `cs_mean`) the returned effect is simply
the posterior coefficient. For covariance-style channels (`corr`,
`vcov`) the returned effects are grouped by their labelled covariance
component.

## Usage

``` r
assoc(object, ...)

posterior_assoc(object, ...)

# S3 method for class 'JoiNMeFit'
assoc(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- ...:

  Unused.

- draws:

  Optional number of posterior draws to retain.

- seed:

  Random seed used when subsetting posterior draws.

- digits:

  Number of digits used when `summary = TRUE`.

- summary:

  Logical. If `TRUE`, return posterior summaries. If `FALSE`, return raw
  MCMC sample matrices.

## Value

A named list with class `PosteriorAssoc`. Each list element contains
either a posterior summary table or a draws-by-term matrix for one
association channel.
