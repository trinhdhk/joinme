# Marker-weight submodel: Transformed parameters

[Back to the marker-weight submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Combines known offsets, common fitted locations and mean-zero unit-scale marker
departures into the weights used by association
channels. For the family name `"student_t"`,
it adds two to one positive Gamma-governed excess to obtain degrees of freedom
shared across all sets, with finite conditional variance.

Effective weights and, where learned, Student-t degrees of freedom are stored
in posterior draws.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
