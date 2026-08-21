# Latent-class submodel: Model

[Back to the latent-class story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Adds priors for baseline probabilities, allocation slopes, and class-specific
locations, then evaluates the selected random-effect blocks under their shared
finite mixture. During dynamic prediction it combines fitted allocation
probabilities with the new subject's longitudinal evidence.

## Files in this block

- [`fit.stan`](fit.stan) — estimation priors and mixture contribution.
- [`dynamic_prediction.stan`](dynamic_prediction.stan) — prediction conditioning contribution.

Reusable probability functions are defined in the preceding
[`functions`](../functions/README.md) stage and imported only by mixture master
programmes.
