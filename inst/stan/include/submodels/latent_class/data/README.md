# Latent-class submodel: Data

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares the number of classes, selected random-effect coordinates, allocation
domains, class-membership design matrices, label-ordering convention, and the
prior information required for baseline probabilities, allocation slopes, and
the within-class component family. The family and any fixed Student-t degrees
of freedom are encoded from `jm_prior(class = list(family = ...))` using the
same family map as the common prior interpreter. Dynamic prediction instead
receives fitted draws and the retained component-family declaration.

## Files in this block

- [`fit.stan`](fit.stan) — estimation data.
- [`dynamic_prediction.stan`](dynamic_prediction.stan) — fitted draws and new-unit allocation designs.
