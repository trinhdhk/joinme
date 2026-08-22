# Latent-class submodel: Functions

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines the finite-mixture component density and class-probability calculation
in separate files. The fitting master imports both definitions. Mixture dynamic
prediction imports only the component density because its class probabilities
arrive as retained posterior draws. Neither ordinary master imports either
definition.

## Files in this block

- [`component_density.stanfunctions`](component_density.stanfunctions) —
  Student-t, Laplace and Normal component densities used in fitting and mixture
  dynamic prediction. Its family code and fixed Student-t degrees of freedom
  come from `jm_priors(class = list(family = ...))`; the fitted class location
  and scale are applied inside this single density function.
- [`class_probability.stanfunctions`](class_probability.stanfunctions) — the
  multinomial-logit class-probability calculation used only in mixture fitting.

The model and generated-quantity fragments call these functions but do not
redeclare them, so the mixture probability law has one source definition.
