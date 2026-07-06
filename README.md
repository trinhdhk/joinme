# Bayesian joint mixed-effects (JoiNMe) model using Stan

`JoiNMe` fits Bayesian joint mixed-effects models for multivariate longitudinal
markers and time-to-event outcomes. It supports multiple outcome families
(Gaussian, Student-t, binary, count, beta, ordinal, skewed families), flexible
association structures, covariance regression, transform-based modelling for
nonlinear effects, and dynamic prediction for individualised trajectory and risk
forecasting.

The package is designed around a single user-facing workflow:

1. simulate or assemble longitudinal and event data,
2. fit the joint model with `joinme()`,
3. inspect renamed posterior draws and diagnostics,
4. summarise and plot fitted effects,
5. generate dynamic predictions and plot future trajectories or survival.

`JoiNMe` works with CmdStanR (preferred) and RStan (compatibility backend). It
also includes helpers for data preparation, posterior diagnostics, simulation,
dynamic prediction, and plotting. See the vignettes and reference for worked
examples and API details.

For marker-aggregated association terms (`cv_total`, `cv_marker`, `cs_total`,
`cs_marker`), transforms are applied at marker level before weighted averaging
in simulation, fitting, and prediction.

Association plotting now follows the fitted model more closely: `plot(fit,
type = "association")` uses a compact draw bundle stored on the fitted object with
the posterior association coefficients, marker weights, spline coefficients,
and cached model-implied raw support ranges. This avoids relying on transient
CmdStan CSV files just to recover association curves and reduces unsupported
support extrapolation for nonlinear transforms.

## Family-shared distributional parameters

When you fit mixed longitudinal families, `JoiNMe` now uses **family-shared**
distributional parameters by default (when no distributional regression is
provided for that parameter).

- Markers with the same family share one latent parameter for each applicable
  distributional quantity.
- Example: for two Student-`t` markers and one Gaussian marker,
  - Student-`t` markers share one `sigma` and one `nu` process,
  - Gaussian marker has its own `sigma` process,
  - `nu` is only defined for Student-`t`.
- This is family-level pooling, not marker-level duplication.

If you specify `formulaDist` for a distributional parameter (for example
`sigma ~ 1 + time`), `JoiNMe` uses that regression structure instead of a
constant family-shared baseline for that parameter. In that case, summaries
report regression terms (fixed/random effects) for the distributional model.

This behaviour is applied consistently across fitting, standata/stancode,
prediction, and summary extraction methods.

## Family-scoped `formulaDist` syntax

You can now scope distributional regressions by family using square brackets on
the LHS:

- Global (applies to all rows where parameter exists):
  - `sigma ~ 1 + time`
- Family-scoped (separate regression block):
  - `sigma[family=student_t] ~ 1 + time`
  - `sigma[family=gaussian] ~ 1 + x1`

Allowed distributional parameters:

- `sigma`
- `nu`
- `phi`
- `alpha` (aliases: `alpha_skew`, `skew`)
- `phi_beta`
- `tau_sde`

Example:

```r
formulaDist <- list(
  sigma[family=student_t] ~ 1 + time + (1 | id),
  sigma[family=gaussian] ~ 1 + x1,
  nu[family=student_t] ~ 1,
  alpha_skew[family=skew_normal] ~ 1 + x1,
  phi[family=negbin2] ~ 1,
  phi_beta[family=beta] ~ 1,
  tau_sde[family=skew_double_exponential] ~ 1
)
```

## Installation

```r
# install.packages("remotes")
# remotes::install_github("trinhdhk/JoiNMe")
```

`JoiNMe` is typically fastest and easiest to run with CmdStanR.

```r
# Preferred: CmdStanR backend
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()

# Compatibility backend: RStan
# install.packages("rstan")
```

### Backend behaviour and threading

- At runtime, `JoiNMe` prefers CmdStanR when CmdStan is installed.
- If the preferred backend is unavailable in the current session, `JoiNMe`
  falls back to the other available Stan backend and warns.
- The runtime now always uses the threaded Stan programs. Setting
  `threads_per_chain = 1` keeps execution serial while reusing the same
  thread-capable Stan code path.
- Install-time RStan precompilation is optional and is attempted only when
  `JoiNMe_COMPILE_RSTAN=1` is set in the installation environment.

### Public fit, prediction, and posterior interfaces

The canonical public fit and prediction classes are now `JoiNMeFit` and
`JoiNMeDynPred`, while `JoiNMeFit` and `PredJoiNMeFit` remain as compatibility
aliases. Objects carry both class labels so existing S3 methods continue to
dispatch without breaking old code.

Posterior draws are now available through a dedicated renamed-draw interface:

- `draws(fit, ...)` and `draws(pred, ...)` return `posterior`-compatible draws
  with user-facing parameter names,
- `as.array(fit)` and `as.array(pred)` return the same renamed posterior arrays,
- `mcmc_plot()` forwards those renamed draws to `bayesplot`,
- `longitudinal_plot()`, `survival_plot()`, `cumhaz_plot()`,
  `association_plot()`, and `diagnostic_plot()` provide explicit plotting entry
  points alongside the main `plot()` methods.

## Quick example

``` r

library(JoiNMe)

set.seed(2026)

# Simulate a small dataset
sim <- simulate_joinme(
  n_id = 10,
  families = list(
    jm_family("student_t"),
    jm_family("student_t", inv_link = ~ inv_logit(x / 2))
  ),
  n_obs_per_marker_per_id = 4,
  times_obs = seq(0, 5, length.out = 8),
  quadrature_nodes = 31,
  seed = 2026,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6)
)

# Note: when estimating marker weights, compare to sim$truth$marker_weights
# because effective marker intensities are used directly from w_raw inside Stan.

formulaLong <- y ~ 1 + time + x1 +
  (1 + time | id) +
  (0 + x1 + (1 + time | id) | marker)

formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  formulaEvent = formulaEvent,
  dataEvent = sim$dataEvent,
  assoc = c("cv_total"),
  families = list(
    jm_family("student_t"),
    jm_family("student_t", inv_link = ~ inv_logit(x / 2))
  ),
  transforms = list(cv_total = list(type = "identity")),
  control = list(
    parallel_chains = 1,
    iter_warmup = 100,
    iter_sampling = 100,
    seed = 2026,
    refresh = 0,
    adapt_delta = 0.9
  )
)

inherits(fit, "JoiNMeFit")
inherits(fit, "JoiNMeFit")

# Summaries
summary(fit)
diagnosis(fit)
fixef(fit)
coef(fit)

# Renamed posterior draws
draws(fit, variables = c("time", "alpha_cv_total"), format = "draws_df")

# Bayesplot-backed posterior display with renamed variables
mcmc_plot(fit, variable = c("time", "alpha_cv_total"), type = "trace")

# Random effects / covariance / combined coefficients (nested by formula block)
re_fit <- ranef(fit)
cf_fit <- coef(fit)
vc_fit <- vcov(fit)
# Example accessors:
# re_fit$formulaLong$id
# cf_fit$formulaLong$id
# re_fit$formulaDist$sigma$allFamilies
# vc_fit$formulaDist$nu$allFamilies

# Dynamic prediction for one subject
ndL <- sim$dataLong[sim$dataLong$id == 1, ]
ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]
time_start <- max(ndL$time)

pred <- posterior_epred(
  fit,
  newdataLong = ndL,
  newdataEvent = ndE,
  time_start = time_start,
  times = seq(time_start, time_start + 1, length.out = 20),
  control = list(n_samples = 50),
  seed = 69
)

# Prediction diagnostics + conditional random-effect summaries
summary(pred)
diagnosis(pred)

# Available only when marker covariance depends on id
if (isTRUE(pred$metadata$marker_corr_depends_on_id)) {
  re_pred <- ranef(pred)
  vc_pred <- vcov(pred)
}

# Plot longitudinal and survival predictions
plot(pred, type = c("longitudinal", "survival"), combined = TRUE)
longitudinal_plot(pred)
survival_plot(pred)

# Alternative longitudinal display: posterior mean change heatmap
plot(
  fit,
  type = "longitudinal",
  longitudinal_style = "heatmap",
  scale = "epred",
  threshold = 0.05
)

# Conditioned plotting, compatible with brms::conditional_effects workflows
cond <- make_conditions(sim$dataEvent, vars = "x1")
plot(
  fit,
  type = "longitudinal",
  longitudinal_style = "heatmap",
  condition = cond,
  scale = "epred"
)

# Posterior summary/extraction helpers
posterior_summary(fit)
posterior_fixef(fit)
posterior_ranef(fit)
posterior_coef(fit)
posterior_assoc(fit, summary = TRUE)

# Explicit fitted-object plot helpers
association_plot(fit)
diagnostic_plot(fit, type = "rhat")

# All coefficient extractors support summary = TRUE/FALSE
fixef(fit, summary = TRUE)
ranef(fit, summary = FALSE)
coef(fit, summary = FALSE)

# For multiple subjects with combined=TRUE: returns one combined plot per subject
# (named list). If combiner packages are unavailable, falls back to the
# standard subject/outcome nested list.
# New marker levels in newdataLong are rejected; prediction marker levels must
# match training levels. Dynamic prediction summaries marginalize latent
# augmentation noise by averaging across dynpred posterior rows per stored draw.

# Concordance uses a dense per-subject survival grid internally for stable
# dynamic prediction at the requested horizon.
```

Trinh Dong, 2026
