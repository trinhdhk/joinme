# Diagnostic summary for a latent-progress mixture

Extends the inherited diagnostics with chain-specific Hamiltonian energy
behaviour, exploratory associations between energy and model blocks, and
a direct assessment of the ordinary-scale/component-scale trade-off.

## Usage

``` r
# S3 method for class 'JoiNMeMixFit'
diagnosis(object, draws = NULL, seed = 1, digits = 3, include_corr = TRUE, ...)
```

## Arguments

- object:

  A JoiNMe object.

- ...:

  Additional arguments passed to class-specific methods.

## Value

A `JoiNMeMix_diagnosis` object inheriting from `JoiNMe_diagnosis`.
