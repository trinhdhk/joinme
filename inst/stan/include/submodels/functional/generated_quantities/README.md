# Functional-transformation submodel: Generated quantities

[Back to the functional-transformation submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Exposes fitted affine coefficients required for interpretation.

Transformed association trajectories are otherwise evaluated in the joint event calculation.

## Files in this block

- [`dynamic_prediction_calculations.stan`](dynamic_prediction_calculations.stan) — dynamic-prediction output calculations.
- [`dynamic_prediction_declarations.stan`](dynamic_prediction_declarations.stan) — dynamic-prediction output declarations.
- [`fit_calculations.stan`](fit_calculations.stan) — fitted-model output calculations.
- [`fit_declarations.stan`](fit_declarations.stan) — fitted-model output declarations.

Read the `fit_*` pair for estimation and the `dynamic_prediction_*` pair for prediction. Declarations and calculations are separate because Stan requires every declaration to precede the first executable statement. An empty calculation fragment is deliberate: its Doxygen header states where the corresponding joint calculation occurs.
