# Fit Joint Nested Mixed-effects (JoiNMe) model via cmdstanr or rstan

Fits the Joint Nested Mixed-effects (JoiNMe) Stan model using cmdstanr
(default) or rstan. The engine can be set globally via
options(stan_preferred_engine = "cmdstanr"\|"rstan") or via
control\$engine.

## Usage

``` r
joinme(
  formulaLong,
  dataLong,
  formulaEvent,
  dataEvent,
  formulaVCov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = joinme_priors(),
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  basehaz = joinme_basehaz(),
  seed = NULL,
  ...
)
```

## Arguments

- formulaLong:

  Longitudinal formula defining fixed effects, id-level effects, and the
  marker block. The marker block may optionally include an inner
  `( ... | id )` term for marker-by-id random effects. When omitted,
  marker-by-id effects are disabled (Q_idm = 0). Grouping terms may use
  `||` at either level: outer `( ... || marker )` keeps marker-only and
  marker-by-id blocks independent, while inner `( ... || id )` keeps the
  marker-by-id covariance diagonal. Grouping terms may also use
  `weighted(group, weights = <column>)` to declare formula-scoped
  subject/group weights.

- dataLong:

  Long-format longitudinal data with columns for id, marker, time,
  outcome, and covariates referenced in `formulaLong`.

- formulaEvent:

  Survival formula for baseline covariates and event model. Supported
  LHS forms are `survival::Surv(time, status)`,
  `survival::Surv(start, stop, status)`,
  `survival::Surv(time, status, type = "left")`, and
  `survival::Surv(time1, time2, type = "interval2")`. The legacy
  `type = "interval"` representation is not supported.

- dataEvent:

  Event-process data with either one row per id (`Surv(time, status)`)
  or multiple interval rows per id (`Surv(start, stop, status)`).
  Covariates in `formulaEvent` may vary by interval.

- formulaVCov:

  Covariance regression formula for id-specific marker-by-id effects. If
  the marker block omits the inner `( ... | id )`, marker-by-id effects
  are absent and covariance-style associations (`corr`, `vcov`) are not
  allowed. When present, `corr` associations use the off-diagonal
  entries of the subject-specific Cholesky-correlation factor `K`,
  whereas `vcov` associations use those same off-diagonal `K` entries
  together with the subject-specific standard deviations. If both `corr`
  and `vcov` are requested, `vcov` is kept and `corr` is ignored with a
  warning. The default `~ 1` remains a supported intercept-only
  covariance regression.

- formulaDist:

  Optional list of formulas for distributional regression. Supported LHS
  parameters are:

  - `sigma` (scale parameter)

  - `nu` (degrees of freedom for Student-t)

  - `phi` (precision for Negative Binomial 2)

  - `alpha` (skewness parameter)

  - `kappa` (positive sample-size parameter for the Beta distribution,
    with shapes \\\mu\kappa\\ and \\(1-\mu)\kappa\\)

  - `tau` (quantile/asymmetry parameter in \\(0,1)\\ for the skew double
    exponential distribution)

  Three input styles are supported:

  1.  Named list with RHS-only formulas, e.g.
      `list(sigma = ~ 1 + time)`.

  2.  Unnamed list with explicit LHS, e.g. `list(sigma ~ 1 + time)`.

  3.  Family-scoped LHS using brackets, e.g.
      `list(sigma[family=student_t] ~ 1 + time, sigma[family=gaussian] ~ 1 + x1)`.

  Semantics of scoping:

  - Without `[family=...]`, the formula applies to all marker rows where
    that parameter is defined.

  - With `[family=...]`, the formula is applied only to rows of that
    family; other rows receive zero contribution from that scoped block.

  - Multiple scoped formulas for the same parameter are allowed and are
    estimated jointly with separate coefficients.

  Random-effects terms with `|` are supported, but nested random-effects
  formulas are not. Grouping factors in random-effects terms may use
  `weighted(group, weights = <column>)`.

- control:

  Named list containing sampling configuration.

  - cmdstanr::model\$sample() arguments (e.g., chains, parallel_chains,
    iter_warmup, iter_sampling, seed, refresh, adapt_delta,
    max_treedepth).

  - engine: "cmdstanr" or "rstan". Defaults to
    options(stan_preferred_engine).

  - threads_per_chain: integer; the threaded Stan program is always
    used. `threads_per_chain = 1` keeps execution serial while
    preserving the thread-capable kernel. For engine = "rstan",
    threading uses options(stan.thread = threads_per_chain).

  - grainsize: integer; reduce_sum grainsize for the threaded kernel.
    Defaults to the full subject count when `threads_per_chain = 1`, and
    to
    `max(1, min(n_cores, ceiling(n_id/(4*threads_per_chain*chains))))`
    otherwise.

  - force_recompile: logical; recompile the Stan model if needed.

  - quadrature_nodes: optional positive integer total node target for
    survival integration. Allowed values are exactly 7/15/31/41/51/61.
    Only the node count is passed to Stan; GK nodes/weights are fixed in
    the Stan code.

  - vcov_diag_link: "softplus" or "exp" for covariance regression
    diagonals.

  - `tau_fixed`: optional fixed `tau` in \\(0,1)\\ for the skew double
    exponential distribution. It cannot be combined with a `tau`
    distributional regression.

- draws:

  Optional number of posterior draws used for summaries (not sampling).

- families:

  Marker-specific family specification (optional). Can be a character
  vector of family names aligned to marker order, or a list of
  `jm_family(...)` entries with per-marker links. Supported
  links/inverse-links: `identity`, `log`, `logit`, `probit`, `exp`.

- transforms:

  Transformation specifications for association terms. Prefer declaring
  them with `joinme_tf(...)`; raw named lists remain supported. Fit-time
  functional transforms may also request a free affine shift through
  `intercept = TRUE` and/or `slope = TRUE`, for example
  `joinme_tf(cv_total = ~ expit(x, intercept = TRUE, slope = TRUE))`,
  which is fitted as `expit(iota_1 + iota_2 * x)`. See details.

- priors:

  Prior declaration. Prefer `joinme_priors(...)`; raw named lists with
  components `beta`, `alpha`, `iota`, and `lkj` remain supported.

- fixed_marker_weights:

  Logical; if TRUE, marker weights are fixed at the supplied base
  values. If FALSE, marker-weight perturbations are estimated using the
  family selected by `shrinkage` (0 = Student-t(6), 1 = Laplace, 2 =
  Normal).

- shared_marker_weights:

  Logical; if TRUE, all weighted marker-based association terms share
  one marker-weight structure. If FALSE, each active weighted
  marker-based association term gets its own marker-weight structure.
  When `marker_weights` is a named list, use names `cv_total`,
  `cs_total`, `cv_marker`, and `cs_marker`.

- basehaz:

  An object of class `joinme_basehaz` created by
  [`joinme_basehaz()`](https://trinhdhk.github.io/joinme/reference/joinme_basehaz.md).
  This controls the baseline hazard parameterisation and spline basis.
  See
  [`?joinme_basehaz`](https://trinhdhk.github.io/joinme/reference/joinme_basehaz.md)
  for details.

- ...:

  Additional args passed to joinme_standata().

## Details

Transformations are provided through a single `transforms` argument.
Each element (cv_total, cs_total, corr, vcov) is defined as a list
specifying a type and type-specific fields (see
[`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md)
for details). For covariance-style associations, `corr` acts on the
off-diagonal entries of the subject-specific Cholesky-correlation factor
`K`, while `vcov` acts on those same off-diagonal `K` entries together
with the subject-specific standard deviations. If both `corr` and `vcov`
are requested, `vcov` is kept and `corr` is ignored with a warning.

Distributional regression can be specified in two equivalent forms:

1.  A named list with RHS-only formulas, e.g.
    `list(sigma = ~ 1 + time)`.

2.  An unnamed list with explicit LHS parameter names, e.g.
    `list(sigma ~ 1 + time)`.

The time axis is internally scaled for numerical stability; reported
coefficients are rescaled within Stan so fixed effects remain
interpretable on the original scale.

Marker weights (see
[`joinme_standata()`](https://trinhdhk.github.io/joinme/reference/joinme_standata.md))
are used to form marker-average summaries for both current value (CV)
and current slope (CS) association components. When
`shared_marker_weights = TRUE`, all weighted marker-based association
terms share one marker-weight structure. When
`shared_marker_weights = FALSE`, each active weighted marker-based
association term (`cv_total`, `cs_total`, `cv_marker`, `cs_marker`) gets
its own marker-weight structure. When `fixed_marker_weights = FALSE`,
signed perturbations are estimated around the supplied base weights. The
effective marker intensities are used directly, without additional
normalisation, to scale the corresponding association contribution.

The returned `JoiNMeFit` object stores a compact association plotting
bundle containing only the posterior quantities needed to draw
association curves (`alpha_*`, marker weights, spline coefficients, and
cached model-implied raw support ranges). This keeps association
plotting usable after serialization without needing the full transient
CmdStan CSV outputs.

For CmdStanR fits, `joinme()` also eagerly imports the CSV-backed fit
contents in memory before returning. This mirrors the loading step used
by `cmdstanr::save_object()` so later
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) calls do not rely on
the original CmdStan CSV files remaining on disk.

Formula-scoped subject weighting is supported in random-effect grouping
terms via `weighted(group, weights = <column>)`. For example:

- `(1 + time | weighted(id, weights = id_w))`

- `(0 + x1 + (1 + time | weighted(id, weights = id_w)) | weighted(marker, weights = marker_w))`

- `sigma ~ 1 + (1 | weighted(id, weights = id_w))`

Weighted grouping affects latent random-effect priors for the referenced
grouping factor and applies subject-level likelihood weighting for the
longitudinal and survival contributions. Weight columns must be finite
and strictly positive.

Transformation for supported association terms can be fed to `joinme()`
in a named list. For example
`list(cv_total = list(type = "functional", expr = ~ log1p(x)), cs_total = list(type = "identity"), corr = list(type = "pwlin", x = c(-2, 0, 2), y = c(0.2, 1, 0.2)), vcov = list(type = "identity"))`

Transformation parameterisation:

- Omitted term, `NULL`, or `list(type = "identity")`: identity transform
  (default).

- `list(type = "functional", expr = ~ log1p(x))`: functional transform;
  `expr` is required and has no additional defaults.

- `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
  monotone I-spline; `knots` and `coeff` are required, `degree` defaults
  to `3`.

- `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
  penalised monotone I-spline. This is intended for simulation but it
  works here too, off-labelly; defaults are `n_knots = 6`, `degree = 3`,
  `lambda = 1` when those values are omitted.

- `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
  penalised monotone I-spline with Stan-estimated coefficients; if
  `knots` are omitted they are derived from `x` quantiles using
  `n_knots = 6` by default.

- `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
  monotone I-spline evaluated on `plogis(x)`; explicit `knots` are
  specified on the expit scale because the spline basis itself lives on
  that bounded `expit(x)` domain. This maybe useful for numerical
  stability but highly experimental and may be useful to put for knots
  at the each (near 0 and 1) if users is concerned about the curve at
  the tails.

- `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
  penalised monotone I-spline on `plogis(x)`. This is intended for
  simulation but it works here too, off-labelly.

- `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
  penalised monotone I-spline on `plogis(x)` with Stan-estimated
  coefficients.

- `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`:
  piecewise-linear transform; both `x` and `y` are required. This is
  intended for simulation but it works here too, off-labelly.

For monotone spline transforms:

- `type = "ispline"`: provide `knots` and `coeff` directly (plus
  optional `degree`); `x`/`y`/`lambda` are not used.

- `type = "ispline_penalised"` (alias: `"ispline_penalized"`): provide
  spline structure through `knots` (or `n_knots`) and smoothness penalty
  `lambda`.

  - If `y` is supplied, `JoiNMe` first fits the monotone spline to
    training pairs `(x, y)` in R and passes fixed coefficients to Stan.
    These plug-in coefficients use the same anchored convention as the
    Stan-estimated path: increasing splines run from `0` to `1`, while
    decreasing splines run from `1` to `0`.

  - If `y` is omitted, Stan estimates the monotone spline coefficients
    directly. Use `direction = "decreasing"` when the monotone transform
    should fall as the raw association feature increases.

- `x`: raw association-feature values used either to define training
  pairs or to help derive knot locations.

- `y`: optional target transformed values at those `x` points. Supplying
  `y` activates the legacy plug-in fit; omitting it activates Stan
  estimation.

- `lambda`: smoothness control (larger = smoother transform).

- `direction`: monotone orientation for penalised spline families.
  Accepted values are `"increasing"` and `"decreasing"`. When `y` is
  supplied, the plug-in fit infers the direction from the training pairs
  if you omit it. When `y` is omitted and Stan estimates the spline
  directly, the default remains `"increasing"`.

- `type = "ispline_expit"` / `"ispline_expit_penalised"`: same semantics
  as the I-spline variants above, except the spline basis is built on
  `plogis(x)`. This keeps the spline input on the bounded interval \$(0,
  1)\$ and is often more numerically stable when the raw association
  feature spans a wide range. My experiments showed that this may help
  with a workaround for the Boundary knots but may cause considerable
  suppress at the tails if mis-specificied. Again, this may not affect
  the predictive value at much because the tail of the association term
  is often scarced. However, if you think the distribution of the latent
  association term is heavy-tailed, you may want to put more knots at
  the tails (e.g knots = c(0.001, 0.01, 0.5, 0.99, 0.999)) to fight the
  impact of the [`plogis()`](https://rdrr.io/r/stats/Logistic.html)
  transformation. In the future, a normalisation of the latent
  association term may be added to the model; may be dividing by sqrt of
  the second moment.

Additional arguments are forwarded to
[`joinme_standata()`](https://trinhdhk.github.io/joinme/reference/joinme_standata.md)
(e.g., `assoc`, `basehaz`, `basehaz_degree`, `n_knots`, `time_var`, and
`shrinkage`). Use lme4-style `||` in `formulaLong` to request diagonal
random-effect covariance.

## Examples

``` r
if (FALSE) { # \dontrun{
formulaDist <- list(
  sigma[family=student_t] ~ 1 + time + (1 | id),
  sigma[family=gaussian] ~ 1 + x1,
  nu[family=student_t] ~ 1,
  alpha[family=skew_normal] ~ 1 + x1,
  phi[family=negbin2] ~ 1,
  kappa[family=beta] ~ 1,
  tau[family=skew_double_exponential] ~ 1
)
} # }
```
