# Plot JoiNMe conditional effects

Plots the process-specific effect tables returned by
`conditional_effects(..., plot = FALSE)`. Each process is delegated to
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)'s
plotting method, retaining its line, interval, surface, point, rug, and
theme conventions.

## Usage

``` r
# S3 method for class 'JoiNMeConditionalEffects'
plot(x, plot = TRUE, ask = FALSE, ...)
```

## Arguments

- x:

  A `JoiNMeConditionalEffects` object.

- plot:

  Logical. If `TRUE`, draw each plot in the active graphics device. If
  `FALSE`, return the plots without drawing them.

- ask:

  Logical. Whether to prompt before drawing a subsequent plot.

- ...:

  Arguments passed to the `brms_conditional_effects` plot method,
  including `points`, `rug`, `stype`, `line_args`, `surface_args`, and
  `theme`.

## Value

Invisibly, a named list by process and effect containing `ggplot`
objects.
