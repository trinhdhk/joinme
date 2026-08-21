# Longitudinal submodel: Data

[Back to the longitudinal submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines the observed multivariate marker histories, response families and inverse links, the four nested fixed/random-effect designs, distributional regressions, and the subject-specific covariance-regression inputs. The subject index is established here because it is the common observational axis for the event submodel.

During dynamic prediction the corresponding file supplies the new subject's observed history, future marker grid and retained longitudinal draws.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
