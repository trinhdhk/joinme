# Survival submodel: Generated quantities

[Back to the survival submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares subject-level event log likelihood, cumulative hazard and survival for
a fitted model, plus conditional survival arrays for dynamic prediction.

Calculations reconstruct the joint association predictor before numerical
integration. Fitted quantities loop over all risk rows for a subject and use
the same right-, exact-, left- or interval-censoring contribution as the model.
For a proper interval \((L,R]\), the reported log likelihood includes both
known survival to \(L\) and failure before \(R\).

The reconstructed cumulative hazard uses the node-specific Cox matrix as well
as the node-specific baseline and association designs. Consequently generated
event quantities use the same time-varying Cox predictor as the fitted
likelihood, including in latent-class models.

## Files in this block

- [`dynamic_prediction_calculations.stan`](dynamic_prediction_calculations.stan) — dynamic-prediction output calculations.
- [`dynamic_prediction_declarations.stan`](dynamic_prediction_declarations.stan) — dynamic-prediction output declarations.
- [`fit_calculations.stan`](fit_calculations.stan) — fitted-model output calculations.
- [`fit_declarations.stan`](fit_declarations.stan) — fitted-model output declarations.

Read the `fit_*` pair for estimation and the `dynamic_prediction_*` pair for prediction. Declarations and calculations are separate because Stan requires every declaration to precede the first executable statement. An empty calculation fragment is deliberate: its Doxygen header states where the corresponding joint calculation occurs.
