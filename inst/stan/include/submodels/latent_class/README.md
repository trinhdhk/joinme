# Latent-class submodel

This directory gives the complete Stan story for the latent-progress model.
One class label is shared by every selected random-effect block, so combining
subject, marker, correlation, or covariance components still defines exactly
`n_classes` classes rather than a Cartesian product of memberships.

The baseline class probabilities and the slopes from `formulaClass` determine
allocation probabilities. Class-specific locations shift only the selected
random-effect coordinates. The ordinary longitudinal, survival, association,
marker-weight, and functional submodels remain responsible for their own
scientific likelihood contributions.

The R declaration
`jm_priors(class = list(baseline_prob = ..., slope = ..., family = ...))`
supplies all three prior parts. `baseline_prob` is the Dirichlet concentration,
`slope` is the coefficient prior for `formulaClass`, and `family` is the
centred, unit-scale distribution used within every class. The family is
Student-t, Normal, or Laplace; class-specific locations and scales remain
separate fitted quantities.

## Read the blocks in order

1. [functions](functions/README.md) — component densities and class
   probabilities, included only by mixture masters.
2. [data](data/README.md) — class count, selected coordinates, allocation
   designs, ordering rule, and prior data.
3. [transformed data](transformed_data/README.md) — parameter-free class
   preparation.
4. [parameters](parameters/README.md) — baseline probabilities,
   class-specific locations, and allocation slopes.
5. [transformed parameters](transformed_parameters/README.md) — identified
   locations and regularised class-regression coefficients.
6. [model](model/README.md) — mixture densities, priors, and dynamic-prediction
   conditioning.
7. [generated quantities](generated_quantities/README.md) — posterior class
   probabilities and class-specific summaries.

[Back to the master submodel map](../README.md).
