# Shared Stan mathematics

The scientific declarations and priors now live under
[`submodels`](../submodels/README.md). This directory retains mathematical
functions and genuinely joint calculations which would be misleading if
assigned to one submodel.

- `functions/` contains reusable likelihood, link, association, covariance,
  quadrature, basis and prior functions.
- `bytecode/` contains the reusable scalar and vector expression evaluator.
- `model/fit_threaded_likelihood.stan` and
  `model/dynpred_threaded_likelihood.stan` invoke the joint reduce-sum
  kernels.
- latent-class declarations live under `../submodels/latent_class/`, so
  ordinary entry points do not carry them.
- shared regression-prior declarations and transformations have descriptive
  stage-specific names; the master entry points include them directly between
  the scientific submodel fragments.

A helper file should contain a mathematical operation reused by several
submodels or a calculation whose joint dependence is itself scientifically
meaningful. New component-specific declarations belong under `../submodels/`.
