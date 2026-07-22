# Posterior random-effect alias for fitted JoiNMe models

Posterior random-effect alias for fitted JoiNMe models

## Usage

``` r
posterior_ranef(object, ...)
```

## Arguments

- object:

  A `JoiNMeFit` object.

- ...:

  Additional arguments forwarded to
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

## Value

The same object returned by `ranef(object, summary = FALSE, ...)`.
