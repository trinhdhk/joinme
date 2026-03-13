# Bayesian joint mixed-effects (JoinME) model using Stan

JoinMe fit Bayesian joint mixed-effects modelling of multivariate longitudinal markers and time-to-event outcomes.
It supports multiple outcome families (Gaussian, Student-t, binary, count),
flexible association structures (e.g. conditional or cumulative associations),
covariance regression, transform-based modelling for functional effects, and
dynamic prediction for individualised risk and trajectory forecasts.

JoinMe work with CmdStanR (recommended) or rstan for
estimation, and includes helpers for data preparation, model diagnostics (ELPD,
LOO, WAIC), simulation utilities and plotting. See the vignettes and reference
for worked examples and API details.

For marker-aggregated association terms (`cv_total`, `cv_marker`, `cs_total`,
`cs_marker`), transforms are applied at marker level before weighted averaging
in simulation, fitting, and prediction.

Association plotting now follows the fitted model more closely: `plot(fit,
type = "association")` uses a compact payload stored on the fitted object with
the posterior association coefficients, marker weights, spline coefficients,
and cached model-implied raw support ranges. This avoids relying on transient
CmdStan CSV files just to recover association curves and reduces unsupported
support extrapolation for nonlinear transforms.

## Family-shared distributional parameters

When you fit mixed longitudinal families, `joinme` now uses **family-shared**
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
`sigma ~ 1 + time`), `joinme` uses that regression structure instead of a
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
# remotes::install_github("trinhdhk/joinme")
```

We can use Rstan or cmdstanr.

```r
# CmdStanR backend
# install.packages("cmdstanr", repos = repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# cmdstanr::install_cmdstan()

# RStan backend
# install.packages("rstan")
```

## Quick example

``` r

library(joinme)

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

# Summaries
summary(fit)
diagnosis(fit)

# Random effects / covariance (nested by formula block)
re_fit <- ranef(fit)
vc_fit <- corr(fit)
# Example accessors:
# re_fit$formulaLong$id
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
  vc_pred <- corr(pred)
}

# Plot longitudinal and survival predictions
plot(pred, which = c("longitudinal", "survival"), combined = TRUE)

# For multiple subjects with combined=TRUE: returns one combined plot per subject
# (named list). If combiner packages are unavailable, falls back to the
# standard subject/outcome nested list.
# New marker levels in newdataLong are rejected; prediction marker levels must
# match training levels. Dynamic prediction summaries marginalize latent
# augmentation noise by averaging across dynpred posterior rows per stored draw.
```

Trinh Dong, 2026
