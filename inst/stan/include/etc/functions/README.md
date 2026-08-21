# Shared mathematical functions

These files define reusable Stan functions for links, basis evaluation,
covariance summaries, cumulative hazards, prior families and the joint fit or
prediction kernels. Latent-class functions live with that submodel in
`../../submodels/latent_class/functions/`. Mixture fitting imports its component
density and class-probability definitions; mixture dynamic prediction imports
only the component density. Ordinary masters import neither. Function-level Doxygen
comments state arguments, return values and statistical roles. The large
partial functions are joint by construction: splitting their local state would
obscure the dependence that defines a joint model. Here “partial” means the
subject subset passed to `reduce_sum`; it does not mean a Cox partial
likelihood. The fitting function evaluates the full event-time probability,
including the fitted baseline hazard and censoring contribution.
