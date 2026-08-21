# Stan model map

This directory contains six master entry points:

- `joinme_fit_threading.stan`: ordinary model fitting;
- `joinme_mix_fit_threading.stan`: latent-class model fitting;
- `joinme_dynpred_threading.stan`: ordinary dynamic prediction;
- `joinme_mix_dynpred_threading.stan`: latent-class dynamic prediction;
- `joinme_fitpred_threading.stan`: parameter-free prediction for subjects
  represented in an ordinary fitted model;
- `joinme_mix_fitpred_threading.stan`: the corresponding latent-class route,
  retaining paired fitted random effects and class probabilities.

Each master file states the Stan blocks in their natural order: functions,
data, transformed data, parameters, transformed parameters, model, and
generated quantities. Within every stage it directly includes each scientific
submodel in dependency order. The scientific declarations and calculations are
under [include/submodels](include/submodels/README.md). Shared mathematical
functions and the few common helper cannot be assigned to one submodel
are under [include/etc](include/etc/README.md).

## How to read a complete fitting model

1. Open the selected master entry point.
2. Read its Stan blocks from `data` through `generated quantities`; the
   submodel includes are written directly in their mathematical order.
3. Follow any scientific include to its submodel README and block README.
4. Finish with `include/etc/functions/joinme_fit_partial.stanfunctions`, which
   is the joint likelihood kernel where the submodels meet.

The dynamic-prediction master files follow the same convention. Generated
quantity declarations appear before calculations directly in the master
block. Calculations which condition future survival on the same newly sampled
longitudinal effects are shared helpers because they genuinely combine
submodels.

The fitted-effect prediction master shares the same generated-quantity
equations, but has no parameters or model block. Its data contain realised
subject, marker, marker-by-subject and covariance-factor draws from the fitted
model. This separation makes clear that no second subject-effect distribution
is estimated when `reuse_fitted_re = TRUE`.

## Dependency direction

The longitudinal submodel establishes the common subject, marker and
random-effect dimensions. Marker weighting, survival, association and
functional transformations depend on those dimensions. The master programme
then invokes shared prior transformations, evaluates the coupled likelihood
and, in mixture entry points, evaluates the latent-class submodel. Dependencies should not
point in the opposite direction.
