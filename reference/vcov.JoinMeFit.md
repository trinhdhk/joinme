# Extract posterior covariance summaries

Returns posterior covariance summaries for fitted JoiNMe models.
Supports selecting specific blocks through `what` (for example,
`what = "id"`).

## Usage

``` r
# S3 method for class 'JoiNMeFit'
vcov(object, what = NULL, draws = NULL, ...)
```

## Arguments

- object:

  A JoiNMe fit object.

- what:

  Optional covariance block selector (`"id"`, `"marker"`). When `NULL`
  (default), returns a nested list for all covariance components in
  `formulaLong` and `formulaDist`.

- draws:

  Number of draws to use for summaries.

- ...:

  Unused.

## Value

If `what` is supplied, a data frame for the requested covariance block.
Otherwise, a nested list with `formulaLong` and `formulaDist` covariance
summaries.
