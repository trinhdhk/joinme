# Marker-weight submodel

This directory describes fixed and fitted weights for marker-aggregated
association features.

For weight set $s$ and marker $d$,

$$
\omega_{sd}=\omega^{(0)}_{sd}+\mu_{\omega,s}+z_{sd},
\qquad z_{sd}\sim F_s(0,1).
$$

The offset is declared as `jm_priors(marker_weights = list(offset = ...))` in R.
It may be a wholly named or wholly unnamed vector, with a named collection of
vectors available for unshared association terms. The fitted common location
borrows information across markers and receives the ordinary coefficient prior
declared by `marker_weights$intercept`. The `marker_weights$family` declaration
governs only the standardised marker coordinate, whose location is zero and
ordinary scale is one. That coordinate is the marker-specific departure; there
is no further fitted scale. The association slope multiplies the completed
weighted feature. Its magnitude is informed by marker-to-marker contrast under
the unit-scale departure law, whilst a non-negative orientation removes the
equivalent simultaneous sign reversal of the slope and every weight. Shared
association terms use one set index; unshared terms use separate sets. A common
location for a marker-only channel can remain weakly informed when the average
centred marker trajectory is close to zero, which makes its explicitly declared
prior scientifically important. Under departures selected by the family name
`"student_t"`, all sets share one fitted degrees-of-freedom value

$$
\nu=2+u,\qquad u\sim\operatorname{Gamma}(2,0.1),
$$

where the Gamma distribution uses shape and rate. Thus every marker in every
set informs the same moving tail parameter. A `prior_student_t()` declaration
is a fixed constant instead. The R families `"constant"` and `"none"` omit all
fitted terms and use the offset exactly. They do not require another Stan prior
family code: R sets every relevant parameter dimension to zero and supplies a
valid unused family code for the common prior routine.

## Read the blocks in order

1. [data](data/README.md) — offsets, sharing maps and departure families.
2. [transformed data](transformed_data/README.md) — prepared sharing
   structure.
3. [parameters](parameters/README.md) — common locations, unit-scale marker coordinates,
   optional Student-t tail parameters, and any regularisation auxiliaries.
4. [transformed parameters](transformed_parameters/README.md) — effective
   marker weights.
5. [model](model/README.md) — departure and common-location priors.
6. [generated quantities](generated_quantities/README.md) — why effective
   weights need not be duplicated.

[Back to the master submodel map](../README.md).
