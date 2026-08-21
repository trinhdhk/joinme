# Association submodel

This directory describes how latent longitudinal features enter the event log
hazard.

Association channels include total, subject-mean and marker current values or
slopes, as well as summaries of subject-specific correlation or covariance
matrices. Forward-shifted designs provide finite-difference slopes on the
original time scale.

For a marker-weighted channel, the association slope (`alpha_*` in Stan)
multiplies the marker average after transformation and weighting. Estimated
weights use a unit-scale departure law, so between-marker contrast informs the
slope magnitude. The slope is oriented non-negatively to remove the equivalent
simultaneous sign reversal of the slope and weights. The fitted common weight
location remains inside every effective weight, and the marker-weight submodel
introduces no additional departure-scale multiplier.

## Read the blocks in order

1. [data](data/README.md) — channel activation, current designs and
   forward-shifted designs.
2. [transformed data](transformed_data/README.md) — parameter-free preparation.
3. [parameters](parameters/README.md) — raw association coefficients.
4. [transformed parameters](transformed_parameters/README.md) — effective
   coefficients and the sign convention.
5. [model](model/README.md) — the single application of association priors.
6. [generated quantities](generated_quantities/README.md) — effective
   coefficients and trajectory ownership.

Association curves are evaluated together with the longitudinal predictor and
hazard. R presentation methods then report them by class when the relevant
association block is class-specific.

[Back to the master submodel map](../README.md).
