# Posterior fixed-effect alias for fitted JoiNMe models

Posterior fixed-effect alias for fitted JoiNMe models

## Usage

``` r
posterior_fixef(object, ...)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- ...:

  Additional arguments forwarded to
  [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html).

## Value

The same object returned by `fixef(object, summary = FALSE, ...)`.
