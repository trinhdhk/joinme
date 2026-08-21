# Latent-class submodel: Generated quantities

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Reports posterior membership probabilities for fitted or new allocation units
and preserves class-specific quantities needed by summaries, plots, and
dynamic prediction.

## Files in this block

- [`fit.stan`](fit.stan) — fitted-unit posterior class probabilities and class summaries.
- [`dynamic_prediction_declarations.stan`](dynamic_prediction_declarations.stan) — prediction output declarations.
- [`dynamic_prediction_calculations.stan`](dynamic_prediction_calculations.stan) — prediction probability calculations.
