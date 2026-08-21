# Functional-transformation submodel: Data

[Back to the functional-transformation submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines transformation modes, bytecode, constants, knots and the fitted affine terms used to reshape longitudinal association features.

Dynamic prediction restores the same transformation specification for each fitted draw.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
