# Plot posterior diagnostics from a JoiNMe fit

Plot posterior diagnostics from a JoiNMe fit

## Usage

``` r
diagnostic_plot(
  object,
  type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean",
    "running_quantile"),
  ...
)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- type:

  Diagnostic plot type.

- ...:

  Additional arguments forwarded to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A `ggplot` object.
