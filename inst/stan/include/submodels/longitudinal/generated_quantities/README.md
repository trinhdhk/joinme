# Longitudinal submodel: Generated quantities

[Back to the longitudinal submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares pointwise marker log likelihoods for fitting, and separates
latent-scale, expected-response and predictive-response arrays for dynamic
prediction. Population coefficients and random-effect scales are reported
directly because their sampled values already use the original study-time
basis.

Fit and prediction calculations are separated from declarations where Stan syntax requires every declaration to appear before any executable statement.

## Files in this block

- [`dynamic_prediction_calculations.stan`](dynamic_prediction_calculations.stan) — dynamic-prediction output calculations.
- [`dynamic_prediction_declarations.stan`](dynamic_prediction_declarations.stan) — dynamic-prediction output declarations.
- [`fit_calculations.stan`](fit_calculations.stan) — fitted-model output calculations.
- [`fit_declarations.stan`](fit_declarations.stan) — fitted-model output declarations.

Read the `fit_*` pair for estimation and the `dynamic_prediction_*` pair for prediction. Declarations and calculations are separate because Stan requires every declaration to precede the first executable statement. An empty calculation fragment is deliberate: its Doxygen header states where the corresponding joint calculation occurs.
