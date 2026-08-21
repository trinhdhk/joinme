# Longitudinal submodel: Transformed parameters

[Back to the longitudinal submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Restores coefficients from their prior-standardised form and constructs
subject, marker and subject-by-marker effects together with subject-specific
covariance factors. Every coefficient and covariance factor is expressed in
the longitudinal formula's original study-time basis; no time-dependent row
is rescaled here.

Dynamic prediction reconstructs these quantities locally for one posterior draw, so it has no persistent transformed parameter here.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
