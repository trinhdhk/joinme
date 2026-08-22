# Bayesian joint mixed-effects (JoiNMe) model using Stan

`JoiNMe` (with package name styled as `joinme` for deliberate ambiguity)
fits Bayesian nested mixed-effects models for multivariate longitudinal
markers, either alone or jointly with time-to-event outcomes. It supports
multiple outcome families (Gaussian, Student-t, binary, count, beta, ordinal,
skewed families) and irregular measurement schedules.

Current implementation experiments with flexible association structures, covariance regression, transform-based modelling for nonlinear effects, and dynamic prediction.

`joinme` works with CmdStanR and RStan. The package is designed with recipes for:

1. simulate or assemble longitudinal and event data with `simulate_joinme()`,
2. fit the joint model with `joinme()`,
3. inspect posterior draws and diagnostics,
4. summarise and plot fitted effects,
5. generate dynamic predictions and plot future trajectories or survival.

## Marker distributional parameters

By default, `joinme` uses **family-shared** distributional parameters where possible, i.e.,

- Markers with the same family share one latent parameter for each applicable
  distributional quantity.
- Example: for two Student-`t` markers and one Gaussian marker,
  - Student-`t` markers share one `sigma` and one `nu` process,
  - Gaussian marker has its own `sigma`.

If `formulaDist` is specified for a distributional parameter (for example
`sigma ~ 1 + time`), `joinme` uses that regression structure instead of a
constant family-shared baseline for that parameter. In that case, summaries
report regression terms (fixed/random effects) for the distributional model.
Details below:

## Family-scoped `formulaDist` syntax

You can scope distributional regressions by family using square brackets on
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
- `kappa`: the positive Beta sample size, giving shape parameters
  `mu * kappa` and `(1 - mu) * kappa`
- `tau`: the skew-double-exponential quantile/asymmetry parameter in `(0, 1)`.
  It may be modelled with `formulaDist`, or fixed for one marker with
  `jm_family("skew_laplace", tau = 0.8)`.

Example:

```r
formulaDist <- list(
  sigma[family=student_t] ~ 1 + time + (1 | id),
  sigma[family=gaussian] ~ 1 + x1,
  nu[family=student_t] ~ 1,
  alpha_skew[family=skew_normal] ~ 1 + x1,
  phi[family=negbin2] ~ 1,
  kappa[family=beta] ~ 1,
  tau[family=skew_double_exponential] ~ 1
)
```

For a fixed skew-Laplace quantile, place the constant in the corresponding
marker family rather than in `control`:

```r
families <- list(
  jm_family("gaussian"),
  jm_family("skew_laplace", tau = 0.8)
)
```

The same `families` declaration can be passed to `simulate_joinme()`. The
simulator uses the marker-specific constant to generate outcomes and records it
in the returned distributional truth.

Here `tau = 0.5` gives the symmetric Laplace distribution. Both endpoints are
excluded: `tau = 0` and `tau = 1` are degenerate limits rather than valid
skew-Laplace distributions.

## Installation

```r
# install.packages("remotes")
# remotes::install_github("trinhdhk/joinme")
```

`joinme` is compatible to run with CmdStanR, but it needs compilation.

```r
# Preferred: CmdStanR backend
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()
# precompile_cmdstanr_models()
```

### Backend and threading

- At runtime, `joinme` prefers CmdStanR when CmdStan is installed.
- If the preferred backend is unavailable in the current session, `joinme`
  falls back to the other available Stan backend and warns.

### Public fit, prediction, and posterior interfaces

Posterior draws are available through a dedicated renamed-draw interface:

- `posterior_draws(fit, ...)` and `posterior_draws(pred, ...)` return `posterior`-compatible draws
- `as.array(fit)` and `as.array(pred)` return the same renamed posterior arrays
- `mcmc_plot()` forwards those draws to `bayesplot`, imitating the behaviour of `brms`
- `longitudinal_plot()`, `survival_plot()`, `cumhaz_plot()`,
  `association_plot()`, and `diagnostic_plot()` provide methods for plotting individual processes.f

<!--
Summary-scale consistency:

- Baseline hazard summaries are reported in `summary(fit)$tables$baseline_hazard`
  on both log and hazard-ratio scales.
- Only the baseline-hazard intercept is adjusted by `-log(tmax)`; non-intercept
  basis coefficients are kept as fitted because time scaling contributes a
  constant log-hazard offset.
- Baseline survival covariates are reported separately in
  `summary(fit)$tables$survival_process` with hazard-ratio columns from
  exponentiated log-scale summaries.
-->

## Quick example

```r

library(joinme)

set.seed(1)

# Simulate a small dataset
sim <- simulate_joinme(
  n_id = 10,
  families = list(
    jm_family("student_t"),
    jm_family("student_t", inv_link = ~ inv_logit(x / 2))
  ),
  times_obs = seq(0, 5, length.out = 8),
  quadrature_nodes = 31,
  seed = 2026,
  assoc = c("cv_total"),
  truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.6)))
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

# Summaries
summary(fit)
diagnosis(fit)
fixef(fit)
coef(fit)

# Renamed posterior draws
posterior_draws(fit, variables = c("time", "alpha_cv_total"), format = "draws_df")

# Bayesplot-backed posterior display with renamed variables
mcmc_plot(fit, variable = c("time", "alpha_cv_total"), type = "trace")

# Conditional effects for either joint-model process
profiles <- make_conditions(sim$dataEvent, vars = "x1")
ce <- conditional_effects(
  fit,
  effects = list(longitudinal = "time", event = "x1"),
  conditions = profiles,
  process = c("longitudinal", "event"),
  plot = FALSE
)
plot(ce, ask = FALSE)

# Average marker-specific expected responses after applying each inverse link
ce_average_marker <- conditional_effects(
  fit,
  effects = "time",
  process = "longitudinal",
  longitudinal_estimand = "marginal_marker",
  plot = FALSE
)

# Paired posterior contrast between two covariate profiles
contrast_profiles <- data.frame(
  time = seq(0, 5, length.out = 40),
  x2 = 0,
  cond__ = paste0("time=", round(seq(0, 5, length.out = 40), 2))
)
cc <- conditional_contrast(
  fit,
  groupA = c(x1 = 1),
  groupB = c(x1 = -1),
  conditions = contrast_profiles,
  process = c("longitudinal", "event"),
  longitudinal_estimand = "marginal_marker",
  plot = FALSE
)
plot(cc, condition_variable = "time", ask = FALSE)

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

# Posterior summary/extraction helpers
posterior_summary(fit)
posterior_assoc(fit, summary = TRUE)

# Conventional coefficient interfaces use posterior_summary() for reporting
# and extract() for their draw-level inputs.
fixef(fit)
ranef(fit)
coef(fit)
extract(fit, what = "fixed_effects")
extract(fit, what = "random_effects")
extract(fit, what = "coefficients")

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
```

Trinh Dong, 2026
