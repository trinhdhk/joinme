# Bayes factor between two JoiNMe fits

Bayes factor between two JoiNMe fits

## Usage

``` r
bayes_factor(fit1, fit2, ...)
```

## Arguments

- fit1:

  First fitted `JoiNMeFit` object.

- fit2:

  Second fitted `JoiNMeFit` object.

- ...:

  Additional arguments passed to
  [`bridgesampling::bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html).

## Value

A list with bridge sampling results and Bayes factor.
