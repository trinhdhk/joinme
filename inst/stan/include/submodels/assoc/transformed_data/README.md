# Association submodel: Transformed data

[Back to the association submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Documents that association indexing and shifted designs are completely determined in R.

No additional parameter-free Stan preparation is required.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
