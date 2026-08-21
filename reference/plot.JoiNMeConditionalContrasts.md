# Plot JoiNMe conditional contrasts

Displays the posterior centre and credible interval for each conditional
contrast. A numeric or date-valued condition may be used as a trajectory
axis; otherwise labelled conditions are shown as point intervals.
Markers and competing event types occupy separate panels.

## Usage

``` r
# S3 method for class 'JoiNMeConditionalContrasts'
plot(
  x,
  plot = TRUE,
  ask = FALSE,
  condition_variable = NULL,
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

  A `JoiNMeConditionalContrasts` object.

- plot:

  Logical. If `TRUE`, draw each requested process.

- ask:

  Logical. Whether to prompt before drawing a subsequent process.

- condition_variable:

  Optional single column from `conditions` to place on the horizontal
  axis. If `NULL`, fitted time is preferred when it varies, followed by
  another varying numeric condition and then `cond__`.

- arrange:

  Arrangement of multiple process plots. `"separate"` retains a named
  list and draws one figure at a time; `"grid"` uses an automatically
  sized grid; `"row"` or `"column"` uses a single row or column; and
  `"design"` follows the layout supplied through `design`.

- ncol, nrow:

  Optional positive whole numbers giving the grid dimensions when
  `arrange = "grid"`.

- design:

  A patchwork design string or area description used when
  `arrange = "design"`. Panels follow the process order in `x`.

- widths, heights:

  Optional positive numeric vectors giving relative column widths and
  row heights in the arranged display.

- guides:

  How legends are treated across an arranged display. One of `"keep"`,
  `"collect"`, or `"auto"`.

- ...:

  Unused and reserved for graphical extensions.

## Value

Invisibly, a named list containing one `ggplot` per process when
`arrange = "separate"`, or one `patchwork` display for every other
arrangement.
