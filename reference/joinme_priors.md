# Declare priors by scientific model component

`jm_prior()` names priors by the part of the statistical model they
govern, rather than by internal coefficient letters. `intercept` and
`slope` are global fallbacks. A component-specific declaration replaces
only the named role, leaving the other role to inherit its global
fallback.

Components with ordinary regression roles are `longitudinal`, `vcov`,
and `functional`. `survival` and `assoc` are slope-only because the
baseline hazard supplies the event-process intercept. `marker_weights`
is special: its `offset` gives the known marker-specific contribution,
`shared` selects one common set or term-specific sets, `intercept`
governs each fitted common weight location, and `family` selects only
the centred, unit-scale distribution of marker-specific departures.
Names appearing on the left-hand side of `formulaDist` (`sigma`, `nu`,
`phi`, `alpha`, `kappa`, and `tau`) may be supplied through `...`, each
with its own intercept and slope priors. Bracket selectors refine a
declaration to a family-scoped coefficient block, for example
`` `sigma[family='student']` = list(slope = prior_normal()) ``. A marker
selector such as `` `sigma[marker='y']` `` is accepted when that marker
uniquely identifies one family-scoped coefficient block. If several
markers share the block, use its family selector because the fitted
coefficient is shared and cannot receive different marker-specific
priors.

## Usage

``` r
joinme_priors(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc = NULL,
  functional = NULL,
  marker = NULL,
  lkj = NULL,
  class = NULL,
  .validate = TRUE
)

jm_priors(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc = NULL,
  functional = NULL,
  marker = NULL,
  lkj = NULL,
  class = NULL,
  .validate = TRUE
)

jm_prior(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc = NULL,
  functional = NULL,
  marker = NULL,
  lkj = NULL,
  class = NULL,
  .validate = TRUE
)

joinme_prior(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc = NULL,
  functional = NULL,
  marker = NULL,
  lkj = NULL,
  class = NULL,
  .validate = TRUE
)
```

## Arguments

- ...:

  Named distributional-parameter prior components, such as
  `sigma = list(intercept = ..., slope = ...)`,
  `` `sigma[family='student']` = list(slope = ...) ``, or
  `` `sigma[marker='y']` = list(intercept = ...) ``.

- intercept, slope:

  Global coefficient priors. These arguments must be created with a
  `prior_*()` function. Fixed numerical coefficients belong to
  [`jm_truth()`](https://trinhdhk.github.io/joinme/reference/joinme_truth.md)
  and are rejected here because a constant is not a prior distribution.

- longitudinal:

  A bare prior declaration in complete longitudinal model-matrix order,
  or a named list containing `intercept` and/or `slope` declarations in
  their separate role-wise orders.

- survival:

  Slope prior for event-model covariates. It may be supplied directly or
  as `list(slope = ...)`.

- baseline:

  Intercept and slope priors for a formula-based baseline hazard.

- vcov:

  A shared covariance-regression declaration, an intercept/slope/latent
  list, or a list with `sd` and `corr` components. The last form permits
  distinct population declarations for marginal scales and off-diagonal
  partial correlations. `intercept` and `slope` govern their regression
  coefficients; `latent` governs the loading of unexplained subject
  variation in the corresponding linear predictor.

- marker_weights:

  A list with `offset`, `intercept`, `family`, and `shared`. The logical
  `shared` value determines whether active weighted association terms
  use one common marker-weight set (`TRUE`, the default) or distinct
  sets (`FALSE`). `offset` is the known contribution to each marker
  weight. It accepts one numeric vector whose entries are either all
  named by marker or all unnamed, or a named collection of such vectors
  for unshared association terms. Partly named vectors are rejected
  because their marker alignment is not statistically defined.
  `intercept` is the ordinary coefficient prior for each fitted common
  marker-weight location and may therefore declare its own location and
  scale; when omitted, it inherits the global `intercept` declaration.
  It must be created with a `prior_*()` function. `family` governs only
  the centred unit-scale departure law for marker-specific deviations.
  Accepted stochastic values are `"student_t"`, `"normal"`, `"laplace"`,
  and `"horseshoe"`. Values `"constant"` and `"none"` are synonymous:
  the offset is then used exactly, with no fitted location or departure.
  The family name `"student_t"` fits one degrees-of-freedom value above
  two per active weight set under a shifted `Gamma(2, 0.1)` prior.
  `prior_*()` objects are not accepted for `family` because
  marker-weight departures retain location zero and ordinary scale one.
  `marker_weights` itself must be a named list; bare prior declarations
  are rejected to avoid confusing the common-location prior with the
  departure family. Fixed common locations belong to
  [`jm_truth()`](https://trinhdhk.github.io/joinme/reference/joinme_truth.md).

- assoc:

  Slope prior for longitudinal–event association coefficients.

- functional:

  Intercept and slope priors for fitted affine shifts inside functional
  association transformations.

- marker:

  A family-only prior declaration for standardised marker-level random
  effects. Location zero and scale one are enforced.

- lkj:

  An object from
  [`prior_lkj()`](https://trinhdhk.github.io/joinme/reference/prior_lkj.md)
  or a positive numeric concentration.

- class:

  A named list with `baseline_prob`, the positive Dirichlet
  concentration for baseline class probabilities, and `slope`, the prior
  declaration for coefficients from `formulaClass`. A scalar
  `baseline_prob` is repeated over classes; a vector may provide one
  concentration per class. The slope remains one prior block because
  class-specific formula lists share a reference-class parameterisation.
  When omitted, `class$slope` inherits the global `slope` declaration.

- .validate:

  Logical; if `TRUE` (default), validate the resulting prior
  declarations immediately.

## Value

An object of class `joinme_priors`.

## Examples

``` r
pri <- jm_prior(
  intercept = prior_student_t(df = 6, scale = 2),
  slope = prior_normal(scale = 1),
  longitudinal = list(slope = prior_normal(scale = 0.5)),
  marker_weights = list(
    offset = c(marker_a = 0, marker_b = 0.25),
    intercept = prior_normal(scale = 0.75),
    family = "student_t",
    shared = TRUE
  ),
  assoc = prior_student_t(df = 4, scale = 1),
  sigma = list(intercept = prior_normal(), slope = prior_normal(scale = 0.5)),
  `sigma[family='student']` = list(slope = prior_student_t(df = 4)),
  lkj = prior_lkj(2),
  class = list(
    baseline_prob = c(2, 2, 2),
    slope = prior_normal(scale = 1)
  )
)
print(pri)
#> Prior specification for Joint Mixed Effects model
#>                          component
#>                   global.intercept
#>                       global.slope
#>             longitudinal.intercept
#>                 longitudinal.slope
#>                     survival.slope
#>                 baseline.intercept
#>                     baseline.slope
#>                        assoc.slope
#>               functional.intercept
#>                   functional.slope
#>                  vcov.sd.intercept
#>                      vcov.sd.slope
#>                     vcov.sd.latent
#>                vcov.corr.intercept
#>                    vcov.corr.slope
#>                   vcov.corr.latent
#>              marker_weights.offset
#>              marker_weights.shared
#>           marker_weights.intercept
#>              marker_weights.family
#>                      marker.family
#>                    sigma.intercept
#>                        sigma.slope
#>  sigma[family=student_t].intercept
#>      sigma[family=student_t].slope
#>                       nu.intercept
#>                           nu.slope
#>                      phi.intercept
#>                          phi.slope
#>                    alpha.intercept
#>                        alpha.slope
#>                    kappa.intercept
#>                        kappa.slope
#>                      tau.intercept
#>                          tau.slope
#>                                lkj
#>                class.baseline_prob
#>                        class.slope
#>                                                             value
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                       normal(mu = 0; scale = 0.5)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 1; df = 4)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                                         normal(mu = 0; scale = 1)
#>                                     marker_a = 0, marker_b = 0.25
#>                                                              TRUE
#>                                      normal(mu = 0; scale = 0.75)
#>  student_t(fitted df: 2 + Gamma(2, 0.1); location = 0; scale = 1)
#>                                         normal(mu = 0; scale = 1)
#>                                         normal(mu = 0; scale = 1)
#>                                       normal(mu = 0; scale = 0.5)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 1; df = 4)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                              student_t(mu = 0; scale = 2; df = 6)
#>                                         normal(mu = 0; scale = 1)
#>                                                      lkj(eta = 2)
#>                                                           2, 2, 2
#>                                         normal(mu = 0; scale = 1)
```
