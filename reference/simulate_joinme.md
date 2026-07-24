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

## Usage

``` r
simulate_joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  formulaVCov = ~1,
  formulaDist = NULL,
  formulaAssoc = NULL,
  transforms = NULL,
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  ...,
  seed = .Random.seed[[1]],
  covariate_formulas = list(x1 ~ rnorm(n_id), x2 ~ rnorm(n_id)),
  marker_weights = NULL,
  shared_marker_weights = TRUE,
  fixed_marker_weights = FALSE,
  shrinkage = 0L,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6),
  beta_long = NULL,
  beta_event = NULL,
  dist_coefs = list(),
  re_params = list(id = list(sd = NULL, corr = NULL), marker = list(sd = NULL, corr =
    NULL), id_marker_cov = list(latent = list(sd = NULL, corr = NULL), alpha = NULL, beta
    = NULL, lambda = NULL, diag_link = "softplus"), dist = list()),
  family_params = list(gaussian = list(sigma = 1), student_t = list(sigma = 1.5, nu = 4),
    bernoulli = list(), binomial = list(trials = 10), poisson = list(), negbin2 =
    list(phi = 2), skew_normal = list(sigma = 1, alpha = 0), double_exponential =
    list(sigma = 1), skew_double_exponential = list(sigma = 1, tau = 0.5), beta =
    list(kappa = 10), cumulative_logit = list(cutpoints = c(-1, 1))),
  h0 = NULL,
  baseline_hazard = list(type = "weibull", shape = 1.4, scale = 6),
  formulaBasehaz = NULL,
  beta_basehaz = NULL,
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
  event_var = "event"
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

- formulaVCov:

  Optional covariance-regression formula for the id-specific
  marker-by-id covariance factor (same role as in
  [`joinme_standata()`](https://trinhdhk.github.io/joinme/reference/joinme_standata.md)).
  This formula is evaluated on event-level covariates (one row per
  subject), must not include random-effect bars `( ... | ... )`, and
  must not include the longitudinal time variable.

  Internally, this formula drives subject-specific standard deviations
  and Cholesky-correlation-factor rows used to build `L_i = SD_i * K_i`;
  see `re_params$id_marker_cov`. The default `~ 1` is supported and
  gives an intercept-only covariance regression.

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
    penalised monotone I-spline on `plogis(x)` in legacy plug-in mode.

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

  In other words, simulation currently uses the legacy plug-in spline
  mode; it does not estimate spline coefficients jointly inside Stan.

- n_id:

  Number of subjects.

- families:

  Marker-specific family names. Use
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  entries to supply custom `link`/`inv_link` expressions. A named probit
  link applies `Phi` as its inverse link. In formula expressions,
  `Phi`/`pnorm` are the standard normal CDF and
  `inv_Phi`/`qnorm`/`probit` are the standard normal quantile.
  Simulation evaluates the same bytecode instructions as fitting and
  prediction.

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

- ...:

  Additional arguments passed to fun.

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

- marker_weights:

  Optional base marker weights used for association aggregation. If
  `shared_marker_weights = TRUE`, supply one numeric vector to be shared
  across all active weighted marker-based association terms. If
  `shared_marker_weights = FALSE`, you may instead supply a named list
  with entries `cv_total`, `cs_total`, `cv_marker`, and `cs_marker`.

- shared_marker_weights:

  Logical. If `TRUE`, all active weighted marker-based association terms
  share one marker-weight structure. If `FALSE`, each active weighted
  marker-based association term uses its own marker-weight structure.

- fixed_marker_weights:

  Logical. This has the same meaning as in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).
  If `TRUE`, `marker_weights` are the effective weights and no latent
  perturbation is drawn. If `FALSE` (the default), `marker_weights` are
  base weights and the effective weights are
  `marker_weights + z_marker_weights`. When base weights are omitted
  they default to zero, exactly as they do in
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)
  when marker weights are estimated.

- shrinkage:

  Integer selecting the distribution of each standardized latent
  marker-weight perturbation when `fixed_marker_weights = FALSE`: `0`
  draws Student-t with 6 degrees of freedom, `1` draws standard Laplace,
  and `2` draws standard Normal. This is the same switch used by the
  Stan priors. It does not change explicitly supplied fixed weights.

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
  present, it overrides `assoc`.

- assoc_coefs:

  Association coefficients for hazard terms. Non-`corr`/`vcov` terms
  accept scalar values. The `corr` term accepts a vector of off-diagonal
  `K` coefficients ordered as `(2,1), (3,1), (3,2), ...` in
  lower-triangular row-major order of the marker-by-id random-effect
  covariance dimension. The `vcov` term accepts the same off-diagonal
  `K` entries followed by the subject-specific standard deviations,
  `(2,1), (3,1), (3,2), ..., sd_1, sd_2, ...`; when `||` is used in the
  marker-by-id random-effects block, only the standard deviation entries
  are used.

  Accepted input forms:

  - named numeric vector, e.g. `c(cv_total = 0.4, cs_mean = -0.2)`,

  - named list, e.g.
    `list(cv_total = 0.4, corr = c(0.2, -0.1), vcov = c(0.4, 0.1, 0.5))`.
    Missing channels default to 0.

- beta_long:

  Fixed-effect coefficients for `formulaLong` fixed part. If NULL,
  coefficients are randomly generated and named by model-matrix columns.

- beta_event:

  Survival baseline-covariate coefficients for non-intercept terms in
  `formulaEvent` RHS. If NULL, coefficients are randomly generated.

- dist_coefs:

  Distributional fixed-effect coefficients for `formulaDist` parameters
  (`sigma`, `nu`, `phi`, `alpha`, `kappa`, `tau`).

  For each parameter, coefficients can be:

  - an unnamed numeric vector (matched by column order),

  - a named numeric vector (matched by model-matrix column names),

  - for family-scoped formulas, a named list with per-scope entries.

  Family-scoped list syntax examples:

  - `dist_coefs = list(sigma = list(default = c("(Intercept)" = -0.3), gaussian = c(...), student_t = c(...)))`

  - alias keys like `"sigma[family='student_t']"` are also recognised
    and mapped to the matching family scope.

- re_params:

  Random-effects simulation controls.

  Structure:

  - `id`: controls id-level random effects from `( ... | id)` in
    `formulaLong`.

  - `marker`: controls marker-level random effects from marker-only
    terms.

  - `id_marker_cov`: controls subject-specific covariance-regression for
    marker-by-id latent effects.

  - `dist`: controls random effects for distributional regressions in
    `formulaDist`.

  For `id` and `marker`, each block is a list:

  - `sd`: scalar or length-K vector of random-effect standard
    deviations,

  - `corr`: KxK correlation matrix.

  `id_marker_cov` fields:

  - `latent`: deprecated compatibility input. If supplied, its implied
    lower-triangular factor is folded into the baseline `alpha`
    intercepts before simulation. Marker-by-id latent seeds are still
    drawn as iid standard normal values. Prefer setting `alpha` directly
    in new code.

  - `alpha`: baseline linear predictors for the subject-specific
    covariance regression entries. Diagonal positions control standard
    deviations; off-diagonal positions control the row-wise
    correlation-factor regression, on the tanh scale.

  - `beta`: covariate effects from `formulaVCov` design matrix
    (systematic subject-to-subject covariance shifts by observed
    covariates),

  - `lambda`: non-negative loading on an iid standard-normal subject
    latent perturbation; if a negative value is supplied, the simulator
    folds the sign into the latent draw so the effective model remains
    unchanged but follows the identified convention used during fitting,

  - `diag_link`: link for the subject-specific standard deviations
    (`"softplus"` or `"exp"`).

  Element-wise covariance-regression form is:
  `eta_{im} = alpha_m + x_i^T beta_m + lambda_m z_{im}`, with
  `lambda_m >= 0` and `z_{im} ~ Normal(0, 1)`. Diagonal entries apply
  `diag_link` to give positive subject-specific standard deviations.
  Off-diagonal entries are mapped through `tanh` and then assembled row
  by row into a valid subject-specific Cholesky-correlation factor
  `K_i`. The covariance factor is reconstructed as `L_i = SD_i * K_i`.
  Marker-by-id latent seeds are sampled as iid standard normal values
  and the final marker-by-id effects are obtained as `b_id = L_i z_id`.

  Dimension rules for `id_marker_cov` entries follow marker-by-id
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

- family_params:

  Family-specific constants used when the corresponding parameter has no
  distributional regression. In particular, specify the Beta
  mean/sample-size model with `beta = list(kappa = ...)` and the skew
  double exponential model with
  `skew_double_exponential = list(sigma = ..., tau = ...)`. `kappa` must
  be positive. `tau` must lie in \\(0,1)\\, with \\0.5\\ giving the
  symmetric double exponential distribution.

- h0:

  Optional baseline hazard function `h0(t)` If supplied, it takes
  precedence over `baseline_hazard`/`formulaBasehaz`.

- baseline_hazard:

  Optional baseline hazard specification. Supported forms:

  - character: one of `"constant"`, `"linear"`, `"piecewise"`,
    `"weibull"`, `"spline"`.

  - named list: `list(type = ..., ...)` with mode-specific parameters.

- formulaBasehaz:

  Optional formula-based baseline hazard model on time, e.g.
  `~ 1 + time + I(time^2)`.

- beta_basehaz:

  Optional coefficients for `formulaBasehaz` (aligned by model-matrix
  column names). If NULL, coefficients are generated.

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

- fun:

  Function to apply.

- x:

  Numeric vector on link scale.

- inv_link_bc:

  Parsed inverse-link bytecode.

## Value

A list containing `dataLong`, `dataEvent`, `truth` (also available as
`true_params`), `marker_info`, `helpers`, and `tmax`. Marker-weight
truth distinguishes base, latent, and effective values. Baseline-hazard
truth is stored in `truth$baseline_hazard`; when a log-linear
representation exists, its resolved coefficients are also in
`truth$stan_fit$bs_gamma_c`. Hazard-scale association coefficients are
stored as `alpha_cv_total`, `alpha_cs_total`, `alpha_cv_mean`,
`alpha_cs_mean`, `alpha_corr`, and `alpha_vcov`, matching the fitted
posterior output names.

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
