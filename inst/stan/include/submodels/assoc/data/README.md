# Association submodel: Data

[Back to the association submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines active association channels and the present/forward-shifted longitudinal designs required for current-value and finite-difference slope terms.

Dynamic prediction restores the same channel definitions for each retained draw.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
