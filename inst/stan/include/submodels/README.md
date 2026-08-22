# Stan submodel map

This is the master reading guide for the statistical Stan sources. The
directory has two axes: the first names a scientific submodel and the second
names a Stan programme block. Consequently, the complete story of one
submodel can be read from top to bottom without searching through a collection
of block-wide helper files.

```text
submodels/
├── README.md
├── longitudinal/
│   ├── README.md
│   ├── data/
│   ├── transformed_data/
│   ├── parameters/
│   ├── transformed_parameters/
│   ├── model/
│   └── generated_quantities/
├── survival/                 # same six-block structure
├── assoc/                    # same six-block structure
├── marker_weight/            # same six-block structure
├── functional/               # same six-block structure
└── latent_class/             # same six-block structure
```

Every one of the 36 block directories contains its own `README.md`. Each
block guide states the statistical role of that stage and identifies every
source file in the directory.

## Scientific reading order

1. [Longitudinal model](longitudinal/README.md): multivariate marker outcomes,
   nested random effects, distributional regression and subject-specific
   covariance regression.
2. [Marker weights](marker_weight/README.md): fixed or fitted weights for
   marker-aggregated association features.
3. [Survival model](survival/README.md): cause-specific baseline hazards,
   censoring and quadrature.
4. [Association model](assoc/README.md): current-value, slope, correlation and
   covariance links between the longitudinal and event processes.
5. [Functional transformations](functional/README.md): fixed or fitted
   transformations applied to association features.
6. [Latent classes](latent_class/README.md): shared class membership,
   class-specific random-effect locations, and allocation regression.

The master Stan files two directories above assemble these six stories
directly. Common prior transformations and coupled likelihoods live under
`../etc/` because they combine several submodels rather than constituting a
scientific submodel.

All formula-derived matrices are formed in R on the original study-time
scale. This applies across the longitudinal, survival, distributional,
covariance, class-membership and formula baseline-hazard components. The Stan
programmes scale only event integration coordinates and limits; they consume
the already evaluated original-time matrices, including a separate Cox design
at every quadrature ordinate.

Within any submodel, read the block directories in Stan order:

```text
data → transformed_data → parameters → transformed_parameters
     → model → generated_quantities
```

## File convention inside a block

- `fit.stan` contributes to ordinary model estimation.
- `dynamic_prediction.stan` contributes to ordinary dynamic prediction.
- generated quantities use declaration and calculation files when necessary,
  because Stan requires all declarations before the first executable
  statement.

An empty Stan fragment is intentional. Its Doxygen header explains why that
submodel has no quantity at the corresponding stage and where the relevant
calculation takes place.

## From submodel to executable programme

The fitting and dynamic-prediction master entry points two directories above
contain the complete ordered include story:

| entry point | estimation | latent classes |
|---|---:|---:|
| `joinme_fit_threading.stan` | yes | no |
| `joinme_mix_fit_threading.stan` | yes | yes |
| `joinme_dynpred_threading.stan` | dynamic prediction | no |
| `joinme_mix_dynpred_threading.stan` | dynamic prediction | yes |
| `joinme_fitpred_threading.stan` | fitted-effect prediction | no |
| `joinme_mix_fitpred_threading.stan` | fitted-effect prediction | yes |

