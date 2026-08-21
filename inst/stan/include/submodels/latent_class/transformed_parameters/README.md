# Latent-class submodel: Transformed parameters

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Constructs identified class locations and forms regularised
`formulaClass` slopes. With slope covariates, only the chosen intercept-like
location coordinate is ordered when intercept ordering is requested.

## Files in this block

- [`fit.stan`](fit.stan) — estimation transformations.
- [`dynamic_prediction.stan`](dynamic_prediction.stan) — documented empty prediction contribution.
