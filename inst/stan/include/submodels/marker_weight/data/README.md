# Marker-weight submodel: Data

[Back to the marker-weight submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines offsets supplied through `priors$marker_weights`, sharing maps, fitted common weight locations,
the distribution family for marker-specific departures, and whether the
Student-t degrees of freedom shared across all weights are explicitly fixed or
learned from the markers.
The R families `"constant"` and `"none"` set the fitted dimensions to zero, so
the declared offsets pass directly to the association calculation.

Dynamic prediction receives the fitted effective weights draw by draw.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
