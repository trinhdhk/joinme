# Plot longitudinal trajectories from JoiNMe objects

Plot longitudinal trajectories from JoiNMe objects

## Usage

``` r
longitudinal_plot(
  object,
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit` or `JoiNMeDynPred` object.

- longitudinal_times:

  Optional numeric vector of fitted trajectory times on the original
  study-time scale. This argument applies to `JoiNMeFit` objects;
  existing `JoiNMeDynPred` objects retain their stored time grid.

- longitudinal_points:

  Number of evenly spaced fitted trajectory times used for a `JoiNMeFit`
  object when `longitudinal_times = NULL`.

- ...:

  Additional arguments forwarded to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A `ggplot` object, a combined plot, or a named list of plots.
