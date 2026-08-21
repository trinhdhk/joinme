# Plot posterior association output from a latent-progress fit

[`posterior_assoc()`](https://trinhdhk.github.io/joinme/reference/assoc.md)
and [`assoc()`](https://trinhdhk.github.io/joinme/reference/assoc.md)
continue to return their familiar coefficient tables or draw matrices.
When that object came from
[`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md),
this method can additionally display the fitted association contribution
along each latent-class trajectory.

## Usage

``` r
# S3 method for class 'PosteriorAssoc'
plot(
  x,
  trajectory = TRUE,
  estimand = c("mean_per_class", "marginal_per_class"),
  ...
)
```

## Arguments

- x:

  A `PosteriorAssoc` object returned from a mixture fit.

- trajectory:

  Logical; request the class-specific trajectory display.

- estimand:

  Either `"mean_per_class"` or `"marginal_per_class"`.

- ...:

  Further arguments passed to
  [`association_plot()`](https://trinhdhk.github.io/joinme/reference/association_plot.md).

## Value

A ggplot or combined plot returned by
[`association_plot()`](https://trinhdhk.github.io/joinme/reference/association_plot.md).
