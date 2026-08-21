# Longitudinal submodel: Parameters

[Back to the longitudinal submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares the population coefficients, nested random effects, covariance-regression quantities, distributional parameters and ordinal cutpoints learned from the longitudinal outcomes.

New-subject latent coordinates are deliberately declared in
`include/etc/parameters/dynamic_prediction_subject_effects.stan`: they are shared
by the longitudinal history and event prediction and therefore do not belong
to one scientific submodel.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
