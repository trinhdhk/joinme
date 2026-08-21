# Marker-weight submodel: Parameters

[Back to the marker-weight submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Declares fitted common locations, direct unit-scale marker departures,
optional Student-t
degrees-of-freedom excess, and any horseshoe auxiliaries when weights are
learned. One excess shared by all active sets is declared only when the family name
`"student_t"` requests moving degrees of freedom; `prior_student_t()`
introduces no such parameter.

No additional departure scale is declared: it would be multiplicatively
confounded with the association slope. Fixed weights introduce no parameter;
dynamic prediction conditions on retained effective weights.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
