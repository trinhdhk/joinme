# Marker-weight submodel: Model

[Back to the marker-weight submodel story](../README.md) · [Master submodel guide](../../README.md)

## Statistical purpose

Assigns the chosen mean-zero unit-scale distribution to marker departures and
any family-specific auxiliaries. If the family name `"student_t"` is used, all
sets share `nu = 2 + u`, with `u ~ Gamma(2, 0.1)` in shape--rate form. If `df` is
supplied through `prior_student_t()`, that value is fixed and no
degrees-of-freedom parameter is sampled.

There is no additional fitted standard deviation. The association slope
(`alpha` in Stan and `assoc$slope` in the prior interface) scales the resulting
weighted marker feature in the log hazard. The `iota` coefficients instead
define an optional affine transformation inside `tf()` and are not a
substitute marker-weight scale.

The fitted common location uses the prior declared by
`marker_weights$intercept`. It is carried through the generic intercept slot of
the shared coefficient-prior programme because it is a location coefficient.
The independently declared `marker_weights$family` applies only to the
mean-zero, unit-scale marker departures and is never used as the prior for the
common location.

## Files in this block

- [`dynamic_prediction.stan`](dynamic_prediction.stan) — ordinary dynamic-prediction contribution.
- [`fit.stan`](fit.stan) — ordinary fitted-model contribution.

Read `fit.stan` for estimation and `dynamic_prediction.stan` for prediction. An empty fragment is deliberate: its Doxygen header states why this submodel has no quantity at that stage and where the corresponding calculation occurs.
