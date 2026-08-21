# Latent-class submodel: Parameters

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares baseline class probabilities, class-specific locations for selected
random-effect coordinates, and slopes governing class membership through
`formulaClass`. Estimation includes these unknowns; dynamic prediction treats
their fitted posterior draws as data.

## Files in this block

- [`fit.stan`](fit.stan) — latent-class fitting parameters.
- [`dynamic_prediction.stan`](dynamic_prediction.stan) — documented empty prediction contribution.
