# Extract fixed effects

Extract fixed effects

## Usage

``` r
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

When `summary = TRUE`, a data.frame of posterior summaries. When
`summary = FALSE`, a draws-by-term matrix.
