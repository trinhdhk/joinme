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
sim <- simulate_joinme_joint_student_t_cvtotal(
  n_id = 10,
  D = 2,
  n_t = 4,
  seed = 2026,
  include_marker_only = TRUE
)

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
  families = rep("student_t", 2),
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

# Dynamic prediction for one subject
ndL <- sim$dataLong[sim$dataLong$id == 1, ]
ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]
Tstart <- max(ndL$time)

pred <- posterior_epred(
  fit,
  newdataLong = ndL,
  newdataEvent = ndE,
  Tstart = Tstart,
  times = seq(Tstart, Tstart + 1, length.out = 20),
  n_samples = 50,
  seed = 69
)

# Plot longitudinal and survival predictions
plot(pred)
```

Trinh Dong, 2026
