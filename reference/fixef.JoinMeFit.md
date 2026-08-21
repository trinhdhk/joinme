# Extract fixed effects

Extract fixed effects

## Usage

``` r
.summarise_fixed_effect_posterior(
  object,
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE,
  ...
)

# S3 method for class 'JoiNMeFit'
fixef(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...)
```

## Arguments

- object:

  A JoiNMe fit object.

- draws:

  Number of draws to use for summaries.

- seed:

  Random seed for subsetting draws.

- digits:

  Number of digits to round summary values.

- summary:

  Logical. If `TRUE`, return posterior summaries with the same
  inferential columns used throughout the package. If `FALSE`, return
  the posterior draw matrix.

- ...:

  Unused.

## Value

When `summary = TRUE`, a data frame containing longitudinal population
coefficients, event-process coefficients, and fitted common
marker-weight locations. The `component`, `event`, and `assoc_term`
columns identify their statistical roles; `assoc_term` is particularly
important when marker-weight sets are not shared. When
`summary = FALSE`, a draws-by-term matrix is returned, with event and
term-specific marker-weight identities included in the column labels.
