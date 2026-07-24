# JoiNMe: Association Construction

## 1 Purpose

This vignette focuses on **association term construction**, including
transformation selection, interpretability, and practical trade-offs.

It also explains how those association terms propagate into the public
summary, draw-extraction, and plotting APIs so the same labels are used
consistently in tables, MCMC plots, and fitted association curves.

Association transforms apply to all supported longitudinal families
(including double_exponential, skew_double_exponential, beta, and
cumulative_logit) and to competing-risk survival models.

### 1.1 Step-by-step: how association summaries are constructed

At the event time and at the default 15 Gauss-Kronrod nodes used for
integration (configurable), the longitudinal linear predictor is
decomposed into a mean component and a marker component. From these
components, `JoiNMe` constructs:

- **CV_total**: total current value (mean + marker components).
- **CV_mean**: current value of the mean component only.
- **CV_marker**: current value of the marker component only.
- **CS_total**: total current slope, approximated by a forward finite
  difference and rescaled to the original time units.
- **CS_mean**: slope of the mean component.
- **CS_marker**: slope of the marker component.
- **CORR**: a subject-specific variance summary of the marker-by-id
  covariance structure, constant across time.
- **VCOV**: a subject-specific lower-triangular summary of the
  marker-by-id Cholesky factor `L`, also constant across time. The
  channel uses the entries of `L` directly rather than the reconstructed
  covariance matrix `L L'`.

Each component can be included or excluded using the `assoc` argument.
If `assoc` does not include a component, the corresponding coefficient
is prior-only and does not enter the likelihood.

#### 1.1.1 Scale recovery and exponentiation in summaries

`JoiNMe` reports association and survival-process coefficients on the
linear predictor scale first. Any internal time normalisation used by
the Stan model is recovered on that linear scale before exponentiation.

For multiplicative interpretations:

- `summary(fit)$tables$baseline_hazard` reports baseline-hazard
  coefficients on the log scale and hazard-ratio scale. Only the
  intercept is adjusted by `-log(tmax)`; non-intercept basis
  coefficients are unchanged.
- `summary(fit)$tables$survival_process` reports baseline covariate
  coefficients and `Hazard.Ratio` interval columns from exponentiated
  log-scale summaries.
- `summary(fit)$tables$assoc` reports association coefficients on their
  native linear scale (including marker-weight terms when active).

The intercept-only baseline adjustment follows the proportional-hazards
identity that internal time scaling contributes a constant offset on the
log-hazard scale. Therefore only the intercept is shifted;
spline/formula basis terms are kept as fitted.

The CORR and VCOV components require marker-by-id random effects. If the
marker block does not include an inner `( ... | id )` term, then
`Q_idm = 0` and these covariance-style associations are not available.

## 2 Marker-weighted CV and CS

When multiple biomarkers are modelled, `JoiNMe` constructs
marker-average association components using weights \omega_d (one per
marker). For marker-aggregated terms (`cv_total`, `cv_marker`,
`cs_total`, `cs_marker`), transforms are applied at the marker level
first and then averaged with weights. This matters whenever transforms
are nonlinear.

For example, with marker-specific total current values CV\_{id}(t) and
transform f\_{cv}:

ext{cv\\total}(t) = \frac{1}{D}\sum\_{d=1}^{D} \omega_d\\
f\_{cv}\\\left(CV\_{id}(t)\right),

which differs from f\_{cv}(\frac{1}{D}\sum_d \omega_d CV\_{id}(t))
unless f\_{cv} is linear.

For slopes, let CS\_{id}(t) be a marker-specific slope feature and
f\_{cs} the configured CS transform. Then

ext{cs\\total}(t) = \frac{1}{D}\sum\_{d=1}^{D} \omega_d\\
f\_{cs}\\\left(CS\_{id}(t)\right),

and analogously for `cs_marker`. This is intentionally not equal to
f\_{cs}(\frac{1}{D}\sum_d \omega_d CS\_{id}(t)) unless f\_{cs} is
linear.

In practice, provide fixed weights with `marker_weights` and
`fixed_marker_weights = TRUE`, or set `fixed_marker_weights = FALSE` to
use signed perturbations around base weights, \omega_d =
\omega_d^{(0)} + z_d. The `shrinkage` switch controls the standardized
latent family in both fitting and simulation: `0` is Student-t_6(0,1),
`1` is Laplace(0,1), and `2` is Normal(0,1). With estimated weights and
no supplied base, \omega^{(0)}=0. The simulator stores the base, latent,
and effective truths in `truth$marker_weights_base`,
`truth$marker_weights_latent`, and `truth$marker_weights`; the last of
these is the effective weight used to generate survival times.

### 2.1 Plotting fitted association curves

Code

``` r

association_plot(fit)
```

### 2.2 Posterior association extraction and covariance-style labels

Use `posterior_assoc(fit, summary = TRUE)` (or
`assoc(fit, summary = TRUE)`) to extract association-specific posterior
summaries without requiring the full model summary. Use
`summary = FALSE` to obtain posterior draws per association term.

### 2.3 Renamed posterior draws for association terms

Association terms can be extracted directly on the reporting scale with
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md). The
canonical posterior variables are `alpha_cv_total`, `alpha_cs_total`,
`alpha_cv_mean`, `alpha_cs_mean`, `alpha_corr[...]`, and
`alpha_vcov[...]`. Each is the coefficient actually used in the event
hazard, which is also the scale reported by
[`summary()`](https://rdrr.io/r/base/summary.html) and used by
prediction and plotting methods.

Code

``` r

# All association terms with reporting-friendly names
assoc_draws <- draws(fit, variables = "^alpha_", regex = TRUE, format = "draws_df")
head(assoc_draws)

# Bayesplot visualisation with the same labels used by summary(fit)
mcmc_plot(fit, variable = "alpha_cv_total", type = "dens_overlay")
```

### 2.4 Conditional effects in the two joint-model processes

Conditional effects are distinct from subject-specific dynamic
prediction.
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
varies named predictors while holding other predictors at explicit
profile values and can return the longitudinal process, the event
process, or both. The longitudinal estimand is the population-level
expected marker response by default. Set
`longitudinal_estimand = "marker"` to retain the fitted marker-level
deviation, or `"marginal_marker"` to average the marker-specific
predictions within each posterior draw. Because marker averaging occurs
after each inverse-link transformation, it targets an average expected
response rather than an expected response at an average marker effect.

The event estimand is the direct relative-hazard component
`exp(W * gamma)`; the baseline hazard and longitudinal association
contribution are held fixed when making event-covariate comparisons.

Code

``` r

profiles <- make_conditions(fit$dataEvent, vars = "treatment")

ce <- conditional_effects(
  fit,
  effects = list(
    longitudinal = "time:treatment",
    event = "treatment"
  ),
  conditions = profiles,
  process = c("longitudinal", "event"),
  longitudinal_estimand = "marker",
  plot = FALSE
)

# The nested tables use the same estimate__, se__, lower__, upper__, and cond__
# columns as brms conditional effects.
head(ce$longitudinal[[1]])
head(ce$event[[1]])
plot(ce, ask = FALSE)
```

## 3 Association coefficients and naming

Internally, association coefficients are denoted with the \alpha prefix
(e.g., `alpha_cv_total`). This follows joint-model conventions and
avoids confusion with the linear predictor \eta used throughout the
longitudinal and survival submodels. Total, mean, and marker
associations use a positive non-centred parameterisation \alpha =
z\_{\alpha} \cdot sd\_{\alpha} for stable sampling.

When an association component is not selected in `assoc`, its
coefficient is prior-only and does not enter the likelihood. Summaries
only report active components to keep interpretation focused on the
fitted association structure.

## 4 Transform modes

`joinme` supports the following transformation families for association
features:

- **Identity** (`type = "identity"`): no transformation.
- **Functional** (`type = "functional"`): a user-defined expression
  converted into a fixed bytecode program.
- **Monotone I-spline** (`type = "ispline"` or `"ispline_penalised"`): a
  smooth monotone mapping.
- **Bounded-domain monotone I-spline** (`type = "ispline_expit"` or
  `"ispline_expit_penalised"`): the same monotone spline machinery, but
  applied to `plogis(x)` rather than directly to `x`.
- **Ordered piecewise linear** (`type = "pwlin"`): posterior
  interpolation across fixed knots with monotone knot ordinates
  estimated in Stan.

Interpretation note for recovery experiments:

- If you use a fixed spline (`type = "ispline"` or an expit variant with
  explicit `coeff`), the association coefficient is interpretable on
  that fixed transformed scale and can be checked directly in
  simulate-fit recovery tests.
- If you use a jointly estimated penalised spline (`*_penalised`), the
  transform shape and the association coefficient are learned together.
  The fitted alpha is still meaningful inside the model, but coefficient
  recovery is typically weaker than in the fixed-transform case because
  the data must learn both the transformed scale and the hazard weight
  simultaneously.

### 4.1 Default values and minimal parameterisations

- **Identity**:
  - Syntax: `list(type = "identity")`
  - Default when a transform term is omitted or set to `NULL`.
  - No other parameters are used.
- **Functional**:
  - Syntax: `list(type = "functional", expr = ~ log1p(x))`
  - Required field: `expr`.
  - No other defaults beyond the supplied expression.
- **I-spline**:
  - Syntax:
    `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3))`
  - Required fields: `knots`, `coeff`.
  - Default: `degree = 3` if omitted.
- **Penalised I-spline**:
  - Syntax:
    `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`
  - Defaults: `n_knots = 6`, `degree = 3`, `lambda = 1`.
  - If `knots` are omitted in Stan-estimated mode, they are derived from
    `x` quantiles.
- **Expit-based penalised I-spline**:
  - Syntax:
    `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`
  - Defaults: `n_knots = 6`, `degree = 3`, `lambda = 1`.
  - The raw association feature is mapped through `plogis(x)` at
    runtime, and the spline basis is built on that bounded domain.
    Supply training `x` and explicit `knots` directly on the expit scale
    in `[0, 1]` so the fitted spline domain matches the runtime basis.
- **Ordered piecewise linear**:
  - Syntax:
    `list(type = "pwlin", knots = c(-2, -1, 0, 1, 2), direction = "increasing")`.
  - Required field: `knots` (aliases: `cutpoints` and legacy `x`).
  - Default: `direction = "increasing"`; use `"decreasing"` for a
    falling transform.
  - An optional non-negative `lambda` penalises second differences of
    the ordered knot ordinates; it defaults to zero.

## 5 Functional Transform

Code

``` r

transforms <- list(
  cv_total = list(
    type = "functional",
    expr = ~ sinh(x) + tanh(x)
  ),
  cs_total = list(
    type = "identity"
  ),
  corr = list(
    type = "identity"
  )
)
```

The functional evaluator supports common arithmetic and functions such
as `log`, `exp`, `sqrt`, `sin`, `cos`, `tanh`, plus aliases
`sigmoid`/`expit`, `softplus`/`log1p_exp`, and helpers like `cbrt` and
`power`.

Code

``` r

transforms <- list(
  cv_total = list(type = "functional", expr = ~ inv_logit(x)),
  cv_mean = list(type = "functional", expr = ~ exp(0.2 * x)),
  cv_marker = list(type = "functional", expr = ~ power(x, 2)),
  cs_total = list(type = "functional", expr = ~ sqrt(abs(x) + 1e-6)),
  corr = list(type = "functional", expr = ~ cbrt(x))
)
```

## 6 Monotone I-Spline

For `type = "ispline"`, we specify the transform shape following
[`splines2::iSpline()`](https://wwenjie.org/splines2/reference/iSpline.html)
convention:

- `knots`: knot locations (first and last are boundary knots),
- `coeff`: non-negative I-spline basis coefficients,
- `degree`: spline degree (usually `3`).

Defaults:

- `degree = 3` if omitted.
- `knots` and `coeff` must be supplied explicitly.

For `type = "ispline_expit"`, the same direct specification applies,
except the I-spline basis is evaluated on `plogis(x)`. The runtime
association feature is still shown on its raw scale, but the transform’s
`knots` must be supplied on the expit scale because that is the spline
domain.

This is highly expremental and meant for stabilising the boundary knots
(with is usually unknown). You may want to standardise the longitudinal
outcome so that most association would not squashed to 0 or 1, as
`plogis(x)` is very flat outside of the range `[-5, 5]`. You may also
want to put more knots as the extreme values of `plogis(x)` to overcome
the flatness of the transform.

Code

``` r

transforms <- list(
  cv_total = list(
    type = "ispline",
    knots = c(-1, 0, 1),
    coeff = c(0, 0.5, 1, 1.2, 1.5, 1.6),
    degree = 3
  ),
  cs_total = list(
    type = "identity"
  ),
  corr = list(
    type = "identity"
  )
)
```

## 7 Penalised Monotone I-Spline

For `type = "ispline_penalised"`, there are now two supported modes:

- **Simulation mode**: provide `x` and `y`, and the package first learns
  a monotone spline in R before passing fixed coefficients into Stan.
  This is meant for the simulation but can be used in model fitting as a
  mean to provide some prior the association curve.
- **Modelling mode**: omit `y`, keep `knots` (or `n_knots`) plus
  `lambda`, and Stan estimates the monotone spline coefficients jointly
  with the rest of the model.

Meaning of the arguments:

- `x`: raw association-feature values, used either as training inputs or
  to help choose knots.
- `y`: optional target transformed values at those `x` points. If
  present, you are telling the package the shape you want to
  approximate.
- `lambda`: smoothness penalty weight (higher = smoother, lower = more
  flexible).

Defaults:

- `n_knots = 6` if `knots` are omitted.
- `degree = 3` if omitted.
- `lambda = 1.0` if omitted.
- In helper-based plug-in fitting, `weights = NULL` means equal weights
  and `diff_order = 2` is the default penalty order.

Code

``` r

transforms <- list(
  cv_total = list(
    type = "identity"
  ),
  cs_total = list(
    type = "identity"
  ),
  corr = list(
    type = "ispline_penalised",
    x = seq(-2, 2, length.out = 40),
    y = exp(seq(-2, 2, length.out = 40)),
    n_knots = 6,
    degree = 3,
    lambda = 1.0
  )
)
```

Code

``` r

transforms <- list(
  corr = list(
    type = "ispline_penalised",
    x = seq(-2, 2, length.out = 40),
    n_knots = 6,
    degree = 3,
    lambda = 1.0
  )
)
```

## 8 Ordered Piecewise-Linear Association

For model fitting, the raw association feature is divided by fixed
knots, but the association values at those knots are not supplied as an
outcome-like `y` vector. With (K) knots, Stan estimates a (K-1) simplex
(). The relative transform ordinate at knot (j) is

f_j = s \sum\_{r=1}^{j-1} \zeta_r, \qquad s = \begin{cases} 1, &
\text{increasing}, \\ -1, & \text{decreasing}. \end{cases}

Thus (f_1=0), (f_K=s), and every adjacent difference has the requested
sign. Linear interpolation is used between knots and the boundary
ordinates are held constant outside the knot range. The survival
contribution is (f(x)): () estimates the total log-hazard span and the
simplex estimates how that span is distributed across intervals. This
follows the monotonic-effect separation of magnitude and simplex shape
described in the [brms monotonic-effects
vignette](https://paulbuerkner.com/brms/articles/brms_monotonic.html)
and implemented by
[`brms::mo()`](https://paulbuerkner.com/brms/reference/mo.html).

The `direction` argument orders the standardised transform (f). As in
the `brms` construction, () remains a signed magnitude parameter, so the
direction of the realised log-hazard contribution (f(x)) also depends on
the posterior sign of ().

The zero first ordinate is an identifiability constraint. A common free
intercept added to every ordinate cannot be distinguished from the log
baseline hazard, so the reported products (f_j) are relative log-hazard
contributions. `summary(fit)$tables$piecewise_ordinates` reports their
posterior contribution and corresponding hazard ratios at every knot.

Code

``` r

transforms <- list(
  cv_total = list(
    type = "identity"
  ),
  cs_total = list(
    type = "identity"
  ),
  corr = list(
    type = "pwlin",
    knots = c(-2, -1, 0, 1, 2),
    direction = "increasing"
  )
)
```

The former fitted declaration `x = ..., y = ...` remains readable for
backwards compatibility. Its `x` values become the knots and `y` is used
only to infer direction when direction is omitted; it does not fix or
train the fitted association. In contrast,
[`simulate_joinme()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme.md)
deliberately keeps the old fixed `x`/`y` interpolation so established
simulation scenarios are unchanged.

## 9 Association Diagram

Code

``` default
graph TD
  CV[CV total] --> FCV[f_cv]
  CS[CS total] --> FCS[f_cs]
  V[Corr] --> FCORR[f_corr]
  FCV --> ETA[eta assoc]
  FCS --> ETA
  FCORR --> ETA
```

Association components mapped into the survival predictor.

## 10 Hierarchical Covariance Regression

Code

``` default
graph TD
  X[Covariates] --> LP[Linear predictor]
  UL[Latent effect] --> LP
  LP --> L[Cholesky factor]
  L --> W[Marker-by-id effects]
  Z[Latent z] --> W
```

Covariance regression structure for marker-by-id effects.

## 11 Fit with Custom Association

Code

``` r

fit <- joinme(
  formulaLong = y ~ 1 + time + x1 +
    (1 + time || id) +
    (0 + x1 + (1 + time || id) || marker),
  dataLong = sim$dataLong,
  formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
  dataEvent = sim$dataEvent,
  transforms = transforms,
  assoc = c("cv_total", "corr", "vcov")
)
```
