# Latent-progress mixtures for nested multivariate joint models

`joinme_mix()` fits the JoiNMe model with a finite mixture on one or
more standardised random-effect blocks. It has a separate entry point
and fitted-object class, whilst deliberately retaining the formula
parser, likelihood, association transforms, prediction programme,
diagnostics, and posterior utilities used by
[`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).

The allocation variables are analytically marginalised in Stan. This is
important because Hamiltonian Monte Carlo cannot sample discrete class
labels. Posterior class probabilities are reconstructed for every
subject or marker in generated quantities.

## Usage

``` r
joinme_mix(
  formulaLong,
  dataLong,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaVCov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = joinme_priors(),
  basehaz = joinme_basehaz(),
  n_classes = 2L,
  formulaClass = ~1,
  class_type = "subject",
  class_dimensions = NULL,
  class_ordering = NULL,
  fit = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- formulaLong, dataLong, formulaVCov, formulaDist, control, draws,
  families, :

  transforms,priors,basehaz,fit, seed Arguments with the same meaning as
  in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).

- formulaEvent:

  Optional survival formula. Supply this together with `dataEvent`, or
  leave both `NULL` for a longitudinal-only mixture.

- dataEvent:

  Optional event-process data.

- n_classes:

  Number of latent classes \\G\\; must be at least two.

- formulaClass:

  A one-sided formula for class-membership covariates, or a named list
  with `subject` and `marker` formulae. A list containing exactly
  `n_classes` formulae gives a separate predictor for every class.
  Covariates must be constant within the corresponding allocation unit.
  The intercept is represented by the baseline class probabilities.

- class_type:

  Character vector naming the random-effect class types. The only
  accepted values are `"subject"`, `"marker"`, `"corr"`, and `"vcov"`.

- class_dimensions:

  Optional named list of integer coordinate indices, one entry per
  selected class type. For a single type, an integer vector is also
  accepted directly.

- class_ordering:

  Identification rule for the common class labels. `"intercept"` orders
  only the first selected random-intercept location; `"probability"`
  orders the baseline class probabilities; `"none"` imposes no order.
  The default is `"intercept"` for a shared class formula and `"none"`
  for a list of `n_classes` class-specific formulae.

- ...:

  Additional arguments passed to
  [`joinme_standata()`](https://trinhdhk.github.io/joinme/reference/joinme_standata.md),
  including `assoc`, `id_var`, `marker_var`, and `time_var`.

## Value

If `fit = TRUE`, a `JoiNMeMixFit` object inheriting from `JoiNMeFit`. If
`fit = FALSE`, a prepared `JoiNMeMixStanData` object.

## Details

The longitudinal linear predictor remains

\$\$ \eta\_{id}(t)=x\_{id}(t)^\top\beta+ z_i(t)^\top u_i+z_d(t)^\top
v_d+ z\_{id}(t)^\top w\_{id}. \$\$

Latent-class modelling changes the distribution of selected
*standardised* random effects, not the structural decomposition above.
For a selected vector \\a_j\\, the ordinary density \\D(a_j;0,1)\\ is
replaced by

\$\$ p(a_j)=\sum\_{g=1}^{G}\pi_g
D\\a_j;\mu_g,\operatorname{diag}(s_g)\\. \$\$

The latent-class component density is declared through
`jm_priors(class = list(family = ...))` using
[`prior_student_t()`](https://trinhdhk.github.io/joinme/reference/prior_student_t.md),
[`prior_normal()`](https://trinhdhk.github.io/joinme/reference/prior_normal.md),
or
[`prior_laplace()`](https://trinhdhk.github.io/joinme/reference/prior_laplace.md).
This choice does not set the ordinary priors for longitudinal,
association, functional, marker-effect, or marker-weight coefficients.
The component locations and scales are estimated. By default, Stan
orders only the first selected random-intercept location. Selected
slopes and other coordinates remain unrestricted. Alternatively,
baseline class probabilities may be ordered, or ordering may be disabled
with `class_ordering`.

The supported class types are:

- `"subject"`: the standardised subject effects underlying \\u_i\\;

- `"marker"`: the standardised marker effects underlying \\v_d\\;

- `"corr"`: the off-diagonal correlation coordinates in the
  marker-by-subject covariance regression;

- `"vcov"`: the lower-triangular variance–covariance coordinates in the
  marker-by-subject covariance regression.

Compatible types share a single class allocation. In particular,
`c("subject", "vcov")` gives one \\G\\-class probability vector per
subject. The construction therefore has \\G\\ classes, rather than the
Cartesian product of separate level-specific classes. If both
subject-indexed and marker-indexed types are selected, they use the same
\\G\\ component labels and mixing proportions but retain the natural
allocation unit of their block.

By default the first two coordinates of each random-effect block define
the two-dimensional progress plane. A one-dimensional block is retained
as a well-defined special case. `class_dimensions` can select other
coordinates. `corr` and `vcov` are alternatives because both act on the
same covariance-regression block.

If `formulaEvent` and `dataEvent` are both `NULL`, a likelihood-neutral
event scaffold is created solely to reuse the common design code. Event
times and event indicators are then set to zero in Stan data, so
inference is based only on the longitudinal process. Association terms
are consequently unavailable in this mode.

## Examples

``` r
if (FALSE) { # \dontrun{
mixed_fit <- joinme_mix(
  y ~ time + (1 + time | id) +
    (1 + time + (1 + time | id) | marker),
  dataLong,
  survival::Surv(time, event) ~ treatment,
  dataEvent,
  n_classes = 3,
  class_type = c("subject", "vcov"),
  formulaClass = ~ treatment + age,
  priors = jm_priors(class = list(
    baseline_prob = rep(2, 3),
    family = prior_student_t(df = 6)
  )),
  assoc = c("cv_total")
)

longitudinal_only_fit <- joinme_mix(
  y ~ time + (1 + time | id) +
    (1 + time + (1 + time | id) | marker),
  dataLong,
  n_classes = 2,
  class_type = "subject"
)
} # }
```
