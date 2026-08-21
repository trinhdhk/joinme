# Association submodel: Transformed parameters

[Back to the association submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Restores the association slopes, fixes inactive channels at zero and applies
the sign convention for weighted marker features. Their slopes are constrained
to be non-negative, choosing one of the two equivalent orientations obtained by
reversing both a slope and its marker weights. The magnitude remains fitted.

Prediction uses the retained effective coefficients directly.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
