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
  fit = TRUE,
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

  Optional survival formula for baseline covariates and the event model.
  Supply it together with `dataEvent`, or leave both `NULL` to fit only
  the nested longitudinal mixed model. Supported LHS forms are
  `survival::Surv(time, status)`, `survival::Surv(start, stop, status)`,
  `survival::Surv(time, status, type = "left")`, and
  `survival::Surv(time1, time2, type = "interval2")`. The former
  `type = "interval"` representation is not supported. An `interval2`
  response must contain one row per id and presently describes one event
  type. Its full likelihood is evaluated as \\S(L)-S(R)\\ for \\L\<T\le
  R\\; the lower inspection limit is not treated as delayed entry. The
  event covariates on that row are held over the represented risk time,
  whilst longitudinal association terms remain time-varying. See the
  model-interpretation vignette for all four censoring contributions.

- dataEvent:

  Optional event-process data with either one row per id
  (`Surv(time, status)`) or multiple interval rows per id
  (`Surv(start, stop, status)`). Covariates in `formulaEvent` may vary
  by interval for a counting-process response. Left- and
  interval-censored responses require one row per id. Leave this and
  `formulaEvent` as `NULL` for longitudinal-only fitting.

- formulaVCov:

  Covariance regression specification for id-specific marker-by-id
  effects. Supply one formula to share its observed covariates, or
  `list(sd = ~ ..., corr = ~ ...)` to model standard deviations and
  off-diagonal partial correlations independently. Intercepts are always
  component-specific parameters and are not duplicated in either design.
  If the marker block omits the inner `( ... | id )`, marker-by-id
  effects are absent and covariance-style associations (`corr`, `vcov`)
  are not allowed. When present, `corr` associations use the
  off-diagonal entries of the subject-specific Cholesky-correlation
  factor `K`, whereas `vcov` associations use those same off-diagonal
  `K` entries together with the subject-specific standard deviations. If
  both `corr` and `vcov` are requested, `vcov` is kept and `corr` is
  ignored with a warning. The default `~ 1` remains a supported
  intercept-only covariance regression.

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

  - vcov_diag_link: "softplus" or "exp" for covariance regression
    diagonals.

- draws:

  Optional number of posterior draws used for summaries (not sampling).

- families:

  Marker-specific family specification (optional). Can be a character
  vector of family names aligned to marker order, or a list of
  `jm_family(...)` entries with per-marker links. A skew-Laplace marker
  may fix its quantile/asymmetry parameter through
  `jm_family("skew_laplace", tau = 0.8)`; otherwise `tau` is estimated.
  Supported named forward links are `identity`, `log`, `logit`,
  `probit`, and `exp`;
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  also accepts an invertible formula link or a directly specified
  inverse-link formula. In formula syntax, `inv_Phi`/`qnorm`/`probit`
  denote the standard normal quantile and `Phi`/`pnorm` denote the
  standard normal CDF. Thus a probit forward link is inverted to `Phi`
  before fitting.

- transforms:

  Transformation specifications for association terms. Prefer declaring
  them with `joinme_tf(...)`; raw named lists remain supported. Fit-time
  functional transforms may also request a free affine shift through
  `intercept = TRUE` and/or `slope = TRUE`, for example
  `joinme_tf(cv_total = ~ expit(x, intercept = TRUE, slope = TRUE))`,
  which is fitted as `expit(iota_1 + iota_2 * x)`. See details.

- priors:

  Prior declaration from
  [`jm_prior()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
  Global `intercept` and `slope` declarations may be replaced
  independently within `longitudinal`, `survival`, `vcov`, `assoc`,
  `functional`, `marker_weights`, and named distributional regressions.
  Distributional names may use family or marker selectors, for example
  `` `sigma[family='student']` `` or `` `sigma[marker='y']` ``; marker
  selection requires a uniquely associated family-scoped `formulaDist`
  block. `marker$family`, `class$slope`, and an LKJ declaration govern
  their distinct structures. The `priors$marker_weights` component
  contains the complete marker-weight declaration. Its `offset` is added
  to fitted weights and must be wholly named or wholly unnamed. Set
  `family = "constant"` (or `"none"`) to use that offset exactly.
  Otherwise the model adds a fitted common location and standardised
  marker-specific departures from the declared family.

- basehaz:

  An object of class `joinme_basehaz` created by
  [`joinme_basehaz()`](https://trinhdhk.github.io/joinme/reference/joinme_basehaz.md).
  This controls the baseline hazard parameterisation and spline basis.
  See
  [`?joinme_basehaz`](https://trinhdhk.github.io/joinme/reference/joinme_basehaz.md)
  for details.

- fit:

  logical; if TRUE, the model is fitted and a `JoiNMeFit` object is
  returned. If FALSE, only the Stan data list is returned.

- seed:

  Optional random seed for reproducibility.

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
`priors$marker_weights$shared = TRUE`, all weighted marker-based
association terms share one marker-weight structure. When
`priors$marker_weights$shared = FALSE`, each active weighted
marker-based association term (`cv_total`, `cs_total`, `cv_marker`,
`cs_marker`) gets its own marker-weight structure. When the family is
stochastic, one common mean per weight set and signed marker-specific
departures are estimated around the offset declared in
`priors$marker_weights$offset`. Each offset vector must be wholly named
or wholly unnamed. The common mean uses
`jm_prior(marker_weights = list(intercept = ...))`; standardised
departures use the centred unit-scale family named by
`marker_weights$family`, for example
`jm_prior(marker_weights = list(family = "laplace"))`. The family
intercept is a fitting prior and does not set a simulation truth. A bare
numeric value retains the ordinary prior-scale shorthand. The
marker-only current-value and slope channels can supply little
information about a common shift when the average centred marker
trajectory is close to zero, so this location prior is substantively
important. The family name `"student_t"` learns set-specific degrees of
freedom; each excess above two has a `Gamma(2, 0.1)` shape–rate prior. A
`prior_student_t(df = ...)` declaration always fixes `df` instead. The
family names `"constant"` and `"none"` instead use the declared offset
as the complete marker weight and fit neither a common mean nor
departures. The effective marker intensities are used directly, without
additional normalisation, to scale the corresponding association
contribution.

The returned `JoiNMeFit` object stores a compact association plotting
bundle containing only the posterior quantities needed to draw
association curves (`alpha_*`, marker weights, spline coefficients, and
cached model-implied raw support ranges). This keeps association
plotting usable after serialization without needing the full transient
CmdStan CSV outputs.

For CmdStanR fits, `joinme()` also eagerly imports the CSV-backed fit
contents in memory before returning. This mirrors the loading operation
used by `cmdstanr::save_object()` so later
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) calls do not rely on
the original CmdStan CSV files remaining on disk.

`joinme()` can also fit the nested multivariate longitudinal model
without an event process. Leave both `formulaEvent` and `dataEvent` as
`NULL`. The longitudinal likelihood and every random-effect block are
retained, whereas the survival likelihood and longitudinal–survival
association are removed. Internally, a likelihood-neutral event scaffold
is used only to reuse the common design and Stan infrastructure. Its
event parameters are not reported as fitted scientific results.

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
`list(cv_total = list(type = "functional", expr = ~ log1p(x)), cs_total = list(type = "identity"), corr = list(type = "pwlin", knots = c(-2, 0, 2), direction = "increasing"), vcov = list(type = "identity"))`

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

- `list(type = "pwlin", knots = c(-2, -1, 0, 1, 2), direction = "increasing")`:
  ordered piecewise-linear association estimated jointly in Stan.

For fitted piecewise-linear associations, the knots divide the raw
association feature into linear intervals. With \$K\$ knots, Stan
estimates a \$K-1\$ simplex. Its cumulative sums give ordered relative
association ordinates running from zero to one for an increasing curve,
or from zero to minus one for a decreasing curve. The signed association
coefficient estimates the total log-hazard span. The first ordinate is
anchored at zero because a free common ordinate would be confounded with
the baseline hazard.

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
  `y` activates the plug-in fit; omitting it activates Stan estimation.

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
longitudinal_fit <- joinme(
  formulaLong = y ~ time + (1 + time | id) +
    (1 + time + (1 + time | id) | marker),
  dataLong = dataLong,
  families = rep("gaussian", 3)
)

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
