# Longitudinal submodel

This directory tells the complete Stan-side story of the multivariate nested
longitudinal model.

For subject $i$, marker $d$, and time $t$, the linear predictor combines
population effects, subject effects, marker effects and subject-by-marker
effects. Marker-specific inverse links and response families define the
observation distribution. Distributional parameters may have separate
regressions. Subject-by-marker effects use a subject-specific covariance
factor whose marginal standard deviations and off-diagonal correlations may
have distinct regressions.

R evaluates the fixed, random-effect, distributional and covariance-regression
model matrices before scaling the event timeline. Thus every spline and
interaction entering this submodel retains the units, knots, boundaries and
contrasts defined by its source data. Stan receives numerical designs and does
not reconstruct formula transformations.

## Read the blocks in order

1. [data](data/README.md) — outcomes, families, designs, dimensions and
   covariance-regression inputs.
2. [transformed data](transformed_data/README.md) — packed covariance
   coordinates fixed before sampling.
3. [parameters](parameters/README.md) — population, nested random-effect,
   distributional and covariance-regression unknowns.
4. [transformed parameters](transformed_parameters/README.md) — effective
   coefficients, realised effects and subject-specific covariance factors on
   the original study-time basis.
5. [model](model/README.md) — support checks and longitudinal prior
   contributions.
6. [generated quantities](generated_quantities/README.md) — effective scales,
   pointwise marker likelihoods and longitudinal predictions.

The marker likelihood is evaluated in the joint threaded calculation because
the same latent effects also form event-process association features. That
placement preserves the joint model; it does not change ownership of the
longitudinal definitions documented here.

## Simulation correspondence

`jm_truth()` separates a data-generating population from a fitting prior. If
an ordinary random-effect block omits `sd`, the simulator draws each
population standard deviation independently from `Exponential(1)`, matching
the `tau_*` statements in `model/fit.stan`. If it omits `corr`, the simulator
draws its Cholesky correlation factor from
`lkj_corr_cholesky(lkj_eta)`, matching `Lcorr_u` and `Lcorr_v`. These
population quantities are drawn once for a simulated data set. The subject,
marker or distributional random effects are subsequently drawn conditionally
at their grouping levels.

[Back to the master submodel map](../README.md).
