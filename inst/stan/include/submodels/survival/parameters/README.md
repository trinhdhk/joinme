# Survival submodel: Parameters

[Back to the survival submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares cause-specific baseline-hazard coefficients and event-covariate slopes.

Prediction conditions on retained fitted values and therefore introduces no new survival-only parameter.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
