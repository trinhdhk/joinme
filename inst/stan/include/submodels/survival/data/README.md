# Survival submodel: Data

[Back to the survival submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Defines cause-specific event covariates, centred baseline-hazard bases,
censoring intervals and Gauss--Kronrod quadrature inputs.

`event_censor_type` distinguishes survival through an interval, an exact event
at its upper limit, failure by an upper limit, and failure within a proper
interval. R supplies a proper interval-censored observation as two risk rows:
survival from zero to its lower inspection limit, followed by failure between
the lower and upper inspection limits. `S_entry` and `S_event` are therefore
the limits of a risk contribution; `S_entry` is not necessarily study entry.

Dynamic prediction supplies the conditioning time and future event-time grids
on the same basis.

The integration bounds remain dimensionless, but every formula transformation
has already been evaluated in original study-time units. In `fit.stan`, `W`
contains endpoint Cox covariates and `W_gk` contains node-specific Cox
covariates. In `dynamic_prediction.stan`, the conditioning and future
`mat_cov_hazard_gk_*` matrices play the same role. Baseline-covariate values
are repeated; time-spline columns are recomputed at each ordinate from the
fitted R model-matrix template.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
