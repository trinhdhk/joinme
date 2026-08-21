# Association submodel: Model

[Back to the association submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Documents that association-slope priors are applied once by the common coefficient-prior programme.

No duplicate density is added in this submodel fragment.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
