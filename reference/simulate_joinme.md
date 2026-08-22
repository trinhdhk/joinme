# Simulate joint model data

Simulator for `JoiNMe` that mirrors the fitting syntax as closely as
possible. The simulator supports:

- multivariate outcomes via marker-level families,

- multi-structure random effects from `formulaLong` (id, marker,
  marker-by-id),

- formula-based covariate generation through user-provided random
  generators,

- association-driven survival via `assoc`/`formulaAssoc` terms
  consistent with the model fit interface.

The implementation is model-matrix based end-to-end, so every simulated
component is generated from the same formula machinery used during
fitting.

Population coefficients are declared through `truth`. A `prior_*()`
object generates one population coefficient vector and a finite numeric
vector fixes that vector exactly. Every realised coefficient is stored
in `truth`.

## Usage

``` r
simulate_joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  formulaVCov = ~1,
  formulaDist = NULL,
  formulaAssoc = NULL,
  transforms = NULL,
  truth = joinme_truth(),
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  seed = .Random.seed[[1]],
  covariate_formulas = list(x1 ~ rnorm(n_id), x2 ~ rnorm(n_id)),
  assoc = c("cv_total"),
  vcov_diag_link = "softplus",
  family_params = list(gaussian = list(sigma = 1), student_t = list(sigma = 1.5, nu = 4),
    bernoulli = list(), binomial = list(trials = 10), poisson = list(), negbin2 =
    list(phi = 2), skew_normal = list(sigma = 1, alpha = 0), double_exponential =
    list(sigma = 1), skew_double_exponential = list(sigma = 1, tau = 0.5), beta =
    list(kappa = 10), cumulative_logit = list(cutpoints = c(-1, 1))),
  time_cens = 8,
  eps_cs = 0.001,
  integration_control = list(rel.tol = 1e-06, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1, t_max = 50, expand = 1.7, max_expand = 60L),
  quadrature_nodes = 15L,
  n_workers = 1L,
  use_mirai = TRUE,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  y_var = "y",
  event_time_var = "time",
  event_var = "event",
  .mixture_specification = NULL
)
```

## Arguments

- formulaLong:

  Longitudinal formula (same role as in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)).
  Grouping terms may use `weighted(group, weights = <column>)` to mirror
  fitting syntax. The referenced weight column should be available in
  generated covariates (for example via `covariate_formulas`). In nested
  marker terms, outer `( ... || marker )` keeps marker-only and
  marker-by-id blocks independent, while inner `( ... || id )` makes the
  marker-by-id covariance diagonal.

- formulaEvent:

  Event/survival formula (same role as in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)).
  Set this to `NULL` to simulate only the nested longitudinal process.
  In that case `dataEvent` is `NULL`, the complete scheduled
  longitudinal history is retained, and `assoc` and `formulaAssoc` must
  be absent.

- formulaVCov:

  Optional covariance-regression specification for the id-specific
  marker-by-id covariance factor. A formula is shared by both covariance
  components; `list(sd = ~ ..., corr = ~ ...)` supplies independent
  observed-covariate regressions with exactly the same syntax as
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)
  and
  [`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md).
  Each formula is evaluated on event-level covariates (one row per
  subject), must not include random-effect bars `( ... | ... )`, and
  must not include the longitudinal time variable.

  The `sd` formula drives subject-specific standard deviations and the
  `corr` formula drives Cholesky-correlation-factor rows used to build
  `L_i = SD_i * K_i`; their population intercepts, slopes and latent
  loadings are declared under `truth$vcov`. The default `~ 1` gives
  intercept-only regressions for both components.

- formulaDist:

  Optional distributional regression formulas (same role as in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)).
  Supported LHS parameters are `sigma`, `nu`, `phi`, `alpha` (aliases:
  `alpha_skew`, `skew`), `kappa`, and `tau`.

  Family-scoped syntax is supported with square brackets:
  `param[family=<name>] ~ ...`, for example
  `sigma[family=student_t] ~ 1 + time`.

  If no bracket is used (e.g. `sigma ~ 1 + time`), the formula applies
  to all rows where that parameter exists. If a bracket is used, it only
  applies to rows of that family and is estimated separately from other
  scopes.

- formulaAssoc:

  Optional formula for association terms in hazard, e.g.
  `~ cv_total + corr + vcov`. If provided, it overrides `assoc`.

- transforms:

  Optional transform specifications for association terms. Supports the
  same structure as `joinme(..., transforms=...)` for `cv_total`,
  `cv_mean`, `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`,
  `vcov`. Simulation applies `cv_total` and `cv_marker` transforms at
  marker level before weighted averaging (aligned with fit/predict Stan
  semantics). Functional transforms support arithmetic and common
  nonlinear functions, including `inv_logit`/`expit`/`sigmoid`, `exp`,
  `log`, `sqrt`, `power`, `cbrt`, `softplus`/`log1p_exp`, trigonometric
  and hyperbolic functions, the standard normal CDF (`Phi`, `pnorm`),
  and the standard normal quantile (`inv_Phi`, `qnorm`, `probit`).

  Quick parameterisation reference:

  - `list(type = "identity")` or omitted term: identity transform
    (default).

  - `list(type = "functional", expr = ~ log1p(x))`: functional
    transform; `expr` is required.

  - `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
    direct monotone I-spline; `degree` defaults to `3` if omitted.

  - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
    penalised monotone I-spline; defaults are `n_knots = 6`,
    `degree = 3`, and `lambda = 1` when omitted.

  - `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
    monotone I-spline evaluated on `plogis(x)`; explicit knots are
    supplied on the expit scale.

  - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
    penalised monotone I-spline on `plogis(x)` in plug-in mode.

  - `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`:
    piecewise-linear transform; `x` and `y` are required. Simulation
    deliberately treats these as fixed interpolation pairs. This differs
    from model fitting, where knot ordinates are estimated as an ordered
    simplex construction and `y` no longer fixes the curve.

  For monotone spline transforms during simulation:

  - `type = "ispline"`: provide `knots` and `coeff` directly (plus
    optional `degree`).

  - `type = "ispline_penalised"` (alias: `"ispline_penalized"`): provide
    training pairs `x` and `y`, plus smoothness penalty `lambda`. The
    simulator fits the monotone spline in R before evaluating the
    transformed association. The fitted plug-in spline uses anchored
    endpoint coefficients matching the Stan-estimated path: increasing
    splines run from `0` to `1`, while decreasing splines run from `1`
    to `0`. If `direction` is omitted, the simulator infers it from the
    supplied `(x, y)` pairs.

  - `type = "ispline_expit"` / `"ispline_expit_penalised"`: same
    semantics as the ordinary I-spline variants above, except the spline
    basis is built on `plogis(x)` for a bounded-domain representation.
    Training `x` values and explicit `knots` for these transform types
    are specified on that bounded expit scale. This is useful when the
    raw association feature has long tails or steep nonlinear effects.

- truth:

  A
  [`jm_truth()`](https://trinhdhk.github.io/joinme/reference/joinme_truth.md)
  declaration for population generation. A numeric intercept or slope
  fixes its population coefficient vector; a `prior_*()` declaration
  draws that vector once. This applies to `longitudinal`, `survival`,
  `baseline`, `assoc_coef`, `vcov`, `functional`, and named
  `formulaDist` components. Its `marker_weights$offset` component
  supplies the marker-specific constant contribution and its
  `marker_weights$family` component generates centred marker-specific
  weight departures, because those departures are otherwise random
  simulation parameters. That field accepts a family name such as
  `"student_t"`, `"normal"`, `"laplace"`, or `"horseshoe"`. The names
  `"constant"` and `"none"` use the offset without a random departure.
  The family name `"student_t"` draws one value shared by all active
  sets as `2 + Gamma(2, 0.1)`, matching fitting. Set
  `marker_weights$family` to `"constant"` or `"none"` when the declared
  offset is the complete marker weight. No common location or
  marker-specific departure is then drawn. Under a stochastic family,
  the effective weight is `offset + marker_weight_mean + departure`,
  where the departure is drawn directly from the declared centred
  unit-scale family. No additional marker-weight scale is used: the
  survival association slope already scales the weighted marker feature.
  or another supported family name.

- n_id:

  Number of subjects.

- families:

  Marker-specific family names or
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  declarations. Use
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  entries to supply custom `link`/`inv_link` expressions or to fix the
  skew-Laplace quantile for a marker, for example
  `jm_family("skew_laplace", tau = 0.8)`. A named probit link applies
  `Phi` as its inverse link. In formula expressions, `Phi`/`pnorm` are
  the standard normal CDF and `inv_Phi`/`qnorm`/`probit` are the
  standard normal quantile. Simulation evaluates the same bytecode
  instructions and fixed-quantile selection as fitting and prediction.

- marker_levels:

  Optional marker names; defaults to `m1`, `m2`, ...

- times_obs:

  Scheduled observation time grid used for every `(id, marker)` before
  optional visit-time jitter and post-event censoring are applied.

- obs_time_noise_sd:

  Optional Gaussian noise SD added independently to each simulated
  observation time after the visit schedule is chosen; noisy visit times
  are clipped to `[0, time_cens]`.

- censor_longitudinal_after_event:

  Logical; if TRUE (default), simulated longitudinal observations are
  truncated at the subject event time. If FALSE, the longitudinal
  schedule may continue after the event time up to `time_cens`.

- left_truncation_max:

  Optional non-negative delayed-entry bound. When \> 0, each subject
  receives a sampled entry in
  `[0, min(left_truncation_max, stop_time))`.

- truncate_longitudinal_before_entry:

  Logical; when TRUE and delayed entry is active
  (`left_truncation_max > 0`), longitudinal rows observed before the
  sampled entry time are removed.

- seed:

  RNG seed.

- covariate_formulas:

  Named or LHS formulas used to generate event-level covariates, e.g.
  `list(x1 ~ rnorm(n_id), x2 ~ rt(n_id, df = 5))`. Time-varying step
  covariates are supported via `time_varyring(fun, n_step, steps, ...)`,
  for example
  `list(x ~ time_varyring(rnorm, c(3, 6), mean = 0, sd = 1))` or
  `list(x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1))`.
  In this mode, each id receives stepwise periods over `[0, time_cens]`,
  and each period value is sampled from `fun`.

- assoc:

  Association components (same names as fit): `cv_total`, `cv_mean`,
  `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.

  Meaning of each channel:

  - `cv_*`: current-value channels from longitudinal trajectories,

  - `cs_*`: current-slope channels from finite differences (`eps_cs`),

  - `corr`: off-diagonal entries of the subject-specific
    Cholesky-correlation factor `K` from marker-by-id random effects.

  - `vcov`: the same off-diagonal `K` entries together with the
    subject-specific standard deviations.

  You can also provide `formulaAssoc = ~ ...` to select channels; when
  present, it overrides `assoc`. Covariance-regression population values
  belong to `truth$vcov`. Its `sd` and `corr` components each accept
  `intercept`, `slope`, and `latent` declarations. Their element-wise
  form is `eta_{im} = alpha_m + x_i^T beta_m + lambda_m z_{im}`, with
  `lambda_m >= 0` and `z_{im} ~ Normal(0, 1)`. Diagonal entries apply
  `diag_link` to give positive subject-specific standard deviations.
  Off-diagonal entries are mapped through `tanh` and then assembled row
  by row into a valid subject-specific Cholesky-correlation factor
  `K_i`. The covariance factor is reconstructed as `L_i = SD_i * K_i`.
  Marker-by-id latent seeds are sampled as iid standard normal values
  and the final marker-by-id effects are obtained as `b_id = L_i z_id`.

  In the ordinary non-mixture model, \\\lambda_m\\ is the conditional
  standard deviation of the unexplained subject heterogeneity on
  covariance-predictor coordinate \\m\\, before applying `diag_link` or
  `tanh`. It is not itself an entry of `L_i`, a covariance, or a
  correlation. With `class_type = "corr"` or `"vcov"`, the selected
  `z_{im}` has a class-specific location and scale. Conditional on class
  \\g\\, its contribution to `eta_{im}` consequently has location
  `lambda_m * mix_location[g, m]` and distributional scale
  `lambda_m * mix_scale[g, m]`.

  The covariance-regression dimension follows the marker-by-id
  random-effect dimension `Q_idm`:

  - if covariance is full: `M = Q_idm * (Q_idm + 1) / 2` lower-tri
    entries,

  - if covariance is forced diagonal (`||` in nested id-marker term):
    `M = Q_idm`.

  `dist` block syntax:

  - `re_params$dist[[param]]` applies to all random-effect terms for
    that distributional parameter,

  - `re_params$dist[[param]]$terms[[j]]` optionally sets term-specific
    controls (same `sd`/`corr` fields as above).

- vcov_diag_link:

  Link applied to covariance-regression scale predictors; either
  `"softplus"` or `"exp"`.

- family_params:

  Family-specific constants used when the corresponding parameter has no
  distributional regression. In particular, specify the Beta
  mean/sample-size model with `beta = list(kappa = ...)` and the skew
  double exponential model with
  `skew_double_exponential = list(sigma = ..., tau = ...)`. `kappa` must
  be positive. `tau` must lie in \\(0,1)\\, with \\0.5\\ giving the
  symmetric double exponential distribution. A marker-specific `tau` in
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  takes precedence over this shared family constant and over a `tau`
  distributional regression for that marker, matching the fitted
  likelihood.

- time_cens:

  Administrative censoring horizon.

- eps_cs:

  Finite-difference step for slope-type associations (`cs_*`).

- integration_control:

  Control list passed to
  [`integrate()`](https://rdrr.io/r/stats/integrate.html).

- root_control:

  Root finding control for inverse-CDF sampling.

- quadrature_nodes:

  Optional Gauss-Kronrod node count for simulation-side survival
  integration. Allowed values: 7, 15, 31, 41, 51, 61.

- n_workers:

  Number of workers to use when `use_mirai = TRUE`.

- use_mirai:

  Logical; when TRUE and `n_workers > 1`, use `mirai` for parallel
  simulation if available. Parallel jobs are seeded deterministically
  from the main `seed` for reproducible simulations.

- id_var, marker_var, time_var, y_var, event_time_var, event_var:

  Column names aligned with
  [`joinme_standata()`](https://trinhdhk.github.io/joinme/reference/joinme_standata.md)
  defaults.

- .mixture_specification:

  Private latent-progress simulation specification assembled by
  [`simulate_joinme_mix()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme_mix.md).
  Users should call
  [`simulate_joinme_mix()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme_mix.md)
  rather than supplying this argument directly.

- fun:

  Function to apply.

- ...:

  Additional arguments passed to fun.

- x:

  Numeric vector on link scale.

- inv_link_bc:

  Parsed inverse-link bytecode.

## Value

A list containing `dataLong`, `dataEvent`, `truth`, `marker_info`,
`helpers`, and `tmax`. Marker-weight truth distinguishes base,
common-location, departure, and effective values. Baseline-hazard truth
is stored in `truth$baseline_hazard`; when a log-linear representation
exists, its resolved coefficients are also in
`truth$stan_fit$bs_gamma_c`. Marker-specific fixed skew-Laplace
quantiles are recorded in `truth$stan_fit$use_tau_fixed` and
`truth$stan_fit$tau_fixed`; realised row-specific quantiles are recorded
in `truth$distributional$rowwise$tau`. Hazard-scale association
coefficients are stored as `alpha_cv_total`, `alpha_cs_total`,
`alpha_cv_mean`, `alpha_cs_mean`, `alpha_corr`, and `alpha_vcov`,
matching the fitted posterior output names. `truth$recovery` contains
the paired fitting entry-point name and a directly reusable argument
list. Thus
`do.call(joinme, c(sim$truth$recovery$arguments, list(seed = 1)))`
recreates the full fitted specification without retyping formulae,
priors, association controls, or covariance settings.

List of results.

Numeric vector on response scale.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_joinme(
  n_id = 30,
  families = list(
    jm_family("gaussian"),
    jm_family("student_t", inv_link = ~ inv_logit(x / 2)),
    jm_family("skew_normal")
  ),
  quadrature_nodes = 31,
  formulaDist = list(
    sigma[family=gaussian] ~ 1 + x1,
    sigma[family=student_t] ~ 1 + time,
    nu[family=student_t] ~ 1,
    alpha_skew[family=skew_normal] ~ 1 + x1
  )
)

sim_split <- simulate_joinme(
  n_id = 20,
  formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
  left_truncation_max = 1.0,
  seed = 99
)

sim_tv <- simulate_joinme(
  n_id = 20,
  formulaEvent = survival::Surv(time_start, time_stop, event) ~ x + x2,
  covariate_formulas = list(
    x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1),
    x2 ~ rnorm(n_id)
  ),
  left_truncation_max = 0.5
)
} # }
```
