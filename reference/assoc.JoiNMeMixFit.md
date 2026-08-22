# Posterior associations from a latent-progress mixture

Posterior associations from a latent-progress mixture

## Usage

``` r
# S3 method for class 'JoiNMeMixFit'
assoc(
  object,
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE,
  trajectory = FALSE,
  estimand = c("mean_per_class", "marginal_per_class"),
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  marginal_samples = 32L,
  class_specific = TRUE,
  class_points = 3L,
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- draws:

  Optional number of posterior draws to retain.

- seed:

  Random seed used when subsetting posterior draws.

- digits:

  Number of digits used when `summary = TRUE`.

- summary:

  Logical. If `TRUE`, return posterior summaries. If `FALSE`, return raw
  MCMC sample matrices.

- trajectory:

  If `TRUE`, return class-specific association trajectories rather than
  coefficient summaries.

- estimand:

  Class trajectory estimand used when `trajectory = TRUE`.

- longitudinal_times, longitudinal_points:

  Time-grid controls.

- marginal_samples:

  Within-component Monte Carlo values per posterior draw for a marginal
  class trajectory.

- class_specific:

  Logical. If `TRUE`, coefficient summaries also carry compact
  class-specific association contributions for terms affected by a
  fitted class-specific random-effect block.

- class_points:

  Number of reference time or covariance-covariate values used in the
  compact class-specific association tables. The plotting method use its
  requested full grid.

- ...:

  Unused.

## Value

A `PosteriorAssoc` object, or a class-trajectory ggplot when
`trajectory = TRUE`. For summary-form mixture output, relevant compact
class contributions are retained in the `class_association` attribute
and printed beneath the common coefficient tables.
