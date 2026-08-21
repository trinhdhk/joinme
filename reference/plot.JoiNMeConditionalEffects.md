# Plot JoiNMe conditional effects

Plots the process-specific effect tables returned by
`conditional_effects(..., plot = FALSE)`. Each process is delegated to
[`brms::conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)'s
plotting method, retaining its line, interval, surface, point, rug, and
theme conventions.

## Usage

``` r
# S3 method for class 'JoiNMeConditionalEffects'
plot(
  x,
  plot = TRUE,
  ask = FALSE,
  arrange = c("separate", "grid", "row", "column", "design"),
  ncol = NULL,
  nrow = NULL,
  design = NULL,
  widths = NULL,
  heights = NULL,
  guides = c("keep", "collect", "auto"),
  ...
)
```

## Arguments

- x:

  A `JoiNMeConditionalEffects` object.

- plot:

  Logical. If `TRUE`, draw each plot in the active graphics device. If
  `FALSE`, return the plots without drawing them.

- ask:

  Logical. Whether to prompt before drawing a subsequent plot.

- arrange:

  Arrangement of multiple plots. `"separate"` retains a named list and
  draws one figure at a time; `"grid"` uses an automatically sized grid;
  `"row"` or `"column"` uses a single row or column; and `"design"`
  follows the layout supplied through `design`.

- ncol, nrow:

  Optional positive whole numbers giving the grid dimensions when
  `arrange = "grid"`.

- design:

  A patchwork design string or area description used when
  `arrange = "design"`. The panels follow their returned order: all
  longitudinal effects first, followed by all event effects.

- widths, heights:

  Optional positive numeric vectors giving relative column widths and
  row heights in the arranged display.

- guides:

  How legends are treated across an arranged display. One of `"keep"`,
  `"collect"`, or `"auto"`.

- ...:

  Arguments passed to the `brms_conditional_effects` plot method,
  including `points`, `rug`, `stype`, `line_args`, `surface_args`, and
  `theme`.

## Value

Invisibly, a named list by process and effect when
`arrange = "separate"`, or one `patchwork` display for every other
arrangement.
