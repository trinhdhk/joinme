# Survival submodel: Model

[Back to the survival submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Contributes priors and curvature penalties for each cause-specific baseline hazard.

Event likelihood and future survival remain joint calculations because their predictors depend on longitudinal effects, marker weights and functional transformations.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
