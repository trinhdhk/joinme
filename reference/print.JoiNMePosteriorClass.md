# Print posterior class probabilities

Print posterior class probabilities

## Usage

``` r
# S3 method for class 'JoiNMePosteriorClass'
print(x, max_rows = 50L, ...)
```

## Arguments

- x:

  A `JoiNMePosteriorClass` object returned by
  [`posterior_class()`](https://trinhdhk.github.io/joinme/reference/posterior_class.md).

- max_rows:

  Maximum number of unit-by-class probability rows printed for each
  allocation domain. The returned object always retains the full table.

- ...:

  Unused.

## Value

Invisibly returns `x`.
