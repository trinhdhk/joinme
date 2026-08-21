# Longitudinal submodel: Model

[Back to the longitudinal submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Checks marker-family support and contributes priors for longitudinal distributional quantities, nested random effects and covariance regression.

The observation likelihood remains in the joint threaded likelihood because the same latent effects enter the event association.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
