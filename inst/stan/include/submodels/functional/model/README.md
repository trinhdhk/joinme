# Functional-transformation submodel: Model

[Back to the functional-transformation submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Contributes shape, simplex, smoothness and auxiliary priors for fitted functional transformations.

Affine coefficient priors are supplied through the common functional intercept/slope roles.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
