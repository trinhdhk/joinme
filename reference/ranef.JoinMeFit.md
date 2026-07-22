# Extract random effects

Extract random effects

## Usage

``` r
# S3 method for class 'JoiNMeFit'
ranef(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...)
```

## Arguments

- object:

  A JoiNMeFit fit object.

- draws:

  Number of draws to use for summaries.

- seed:

  Random seed for subsetting draws.

- digits:

  Number of digits to round summary values.

- summary:

  Logical. If `TRUE`, return posterior summaries. If `FALSE`, return the
  posterior extraction on the coefficient scale.

- ...:

  Unused.

## Value

A nested list with top-level entries `formulaLong` and `formulaDist`.
`formulaLong` contains random-effect summaries for longitudinal model
components (`id`, `marker`, `marker_by_id_latent`, when present).
`formulaDist` contains distributional random-effect summaries organised
by parameter and family scope (e.g., `sigma$student_t`,
`nu$allFamilies`).
