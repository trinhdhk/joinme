# Survival submodel: Transformed data

[Back to the survival submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Records that the event bases, interval indices and quadrature designs arrive fully prepared from R.

The intentionally empty fit and prediction fragments make this absence explicit rather than leaving the stage undocumented.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
