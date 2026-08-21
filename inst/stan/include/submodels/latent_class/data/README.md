# Latent-class submodel: Data

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares the number of classes, selected random-effect coordinates, allocation
domains, class-membership design matrices, label-ordering convention, and the
prior information required for baseline probabilities and allocation slopes.
Dynamic prediction instead receives fitted draws of these quantities.

## Files in this block

- [`fit.stan`](fit.stan) — estimation data.
- [`dynamic_prediction.stan`](dynamic_prediction.stan) — fitted draws and new-unit allocation designs.
