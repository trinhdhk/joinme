# A Example Workflow for Joint Modelling with mutlivariate Longitudinal and Survival Data using joinme

## 1 Introduction

This vignette outlines a workflow for using `joinme` in a research
context, suitable for inclusion in high-impact statistical or clinical
journals.

In this document, we will be 1. Generating synthetic data with known
ground truth. 2. Fitting the joint model. 3. Assessing convergence and
model fit with renamed posterior draws. 4. Performing dynamic
predictions for individual subjects. 5. Reporting results with
consistent summaries and plotting helpers.

### 1.1 Program flow in practice

The `joinme` package follows the below data-processing and model-fitting
sequence:

1.  **Data assembly (R)**: build design matrices, scale time to
    $`[0, 1]`$, and prepare Gauss-Kronrod nodes (default 15,
    configurable) for survival integration.
2.  **Parameter declarations (Stan)**: define fixed effects, random
    effects, distributional regression coefficients, baseline hazard
    parameters, and association coefficients.
3.  **Transformed parameters (Stan)**: rescale time-related coefficients
    and construct the random effects used in the likelihood.
4.  **Model block (Stan)**: apply priors and evaluate the likelihood
    using `reduce_sum`, which splits subjects across threads when
    requested.
5.  **Generated quantities (Stan)**: compute per-observation
    log-likelihoods, fitted values, and prediction summaries for
    posterior diagnostics.
6.  **Post-processing (R)**: expose the posterior through
    [`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md),
    summary helpers,
    [`mcmc_plot()`](https://rdrr.io/pkg/joinme/man/mcmc_plot.html), and
    explicit plotting wrappers such as
    [`longitudinal_plot()`](https://rdrr.io/pkg/joinme/man/longitudinal_plot.html)
    and
    [`survival_plot()`](https://rdrr.io/pkg/joinme/man/plot_helpers.html).

This sequence is important because it clarifies where each modelling
choice enters the inference: formula parsing happens in R, while all
probability calculations occur in Stan.

## 2 Simulation and Data Exploration

Before analysing real data, it is highly recommended to perform a
**parameter recovery simulation**. This confirms that the model is
identifiable given your sample size and sampling frequency.

We use the built-in generalised simulator `simulate_joinme` to generate
a dataset where the hazard depends on the “Total Current Value” of the
biomarkers.

Code

``` r

set.seed(2025)

# Simulate a dataset: 
# - 100 subjects
# - 3 biomarkers
# - ~10 observations per marker per subject
sim <- simulate_joinme(
  n_id = 40,
  families = list(
    jm_family("student_t"),
    jm_family("student_t"),
    jm_family("student_t", inv_link = ~ inv_logit(x / 2))
  ),
  n_obs_per_marker_per_id = 6,
  times_obs = seq(0, 8, length.out = 12),
  quadrature_nodes = 31,
  beta_long = c("(Intercept)" = 1.0, "time" = 0.5, "x1" = 0.4),
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.5),
  fixed_marker_weights = FALSE,
  shrinkage = 0L
)

# Gauss-Kronrod nodes/weights are fixed in Stan; `quadrature_nodes` selects
# the node count used to build the design matrices in R.

# `truth$marker_weights` is the effective base + latent value used in the hazard.
# The decomposition is available as marker_weights_base/marker_weights_latent.

dataLong <- sim$dataLong
dataEvent <- sim$dataEvent

# Optional: generate counting-process event rows with delayed entry.
sim_split <- simulate_joinme(
  n_id = 20,
  formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
  covariate_formulas = list(
    x1 ~ time_varyring(rnorm, c(2, 5), mean = 0, sd = 1),
    x2 ~ rnorm(n_id)
  ),
  left_truncation_max = 1.0,
  seed = 2026
)
dataEvent_split <- sim_split$dataEvent

# Preview the data
head(dataLong)
#>   id marker      time       x1        x2           y
#> 1  1     m1 0.0000000 1.370958 0.2059986 -1.51724626
#> 2  1     m1 0.7272727 1.370958 0.2059986  2.70925844
#> 3  1     m1 1.4545455 1.370958 0.2059986 -1.75833339
#> 4  1     m1 2.1818182 1.370958 0.2059986  0.02561363
#> 5  1     m1 2.9090909 1.370958 0.2059986 -1.71592539
#> 6  1     m1 3.6363636 1.370958 0.2059986 -4.28562285
head(dataEvent)
#>   id         x1         x2     time event time_start time_stop
#> 1  1  1.3709584  0.2059986 8.000000     0          0  8.000000
#> 2  2 -0.5646982 -0.3610573 1.052114     1          0  1.052114
#> 3  3  0.3631284  0.7581632 3.641210     1          0  3.641210
#> 4  4  0.6328626 -0.7267048 6.699508     1          0  6.699508
#> 5  5  0.4042683 -1.3682810 5.019946     1          0  5.019946
#> 6  6 -0.1061245  0.4328180 1.797717     1          0  1.797717
```

Visualise your longitudinal trajectories and survival distribution
before modelling.

Code

``` r

# Visualise trajectories for a subset of IDs
subset_ids <- base::sample(unique(dataLong$id), 6)

ggplot(dataLong %>% filter(id %in% subset_ids), aes(x = time, y = y, color = marker, group = marker)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~id) + 
  theme_bw() +
  labs(title = "Longitudinal Trajectories (Subset)", x = "Time", y = "Biomarker Value")
```

![](joinme-workflow_files/figure-html/explore-plot-1.png)

Code

``` r


# Check event rate
table(dataEvent$event)
#> 
#>  0  1 
#> 10 30
prop.table(table(dataEvent$event))
#> 
#>    0    1 
#> 0.25 0.75
```

Ensure the event rate is sufficient (typically \>10-20 events per
parameter in the survival submodel) to support the complexity of the
proposed association structure.

## 3 Model Specification

We specify a joint model where: 1. **Longitudinal Submodel**: Linear
growth over time with random intercepts and slopes for each
subject-marker combination. 2. **Survival Submodel**: The hazard depends
on baseline covariates (`x1`, `x2`) and the current value (`cv_total`)
of the biomarkers.

#### 3.0.1 Families and distributional regression

Each marker can follow a distinct family (for example Gaussian,
Student-$`t`$, Poisson, negative binomial, Bernoulli, beta, or ordinal).
When a family includes additional parameters (for example $`\sigma`$,
$`\nu`$, or $`\phi`$), those parameters can be modelled through
`formulaDist`. This allows heteroscedasticity or covariate-dependent
dispersion while keeping the mean structure aligned across markers.

#### 3.0.2 Association transformations

The association terms can be transformed using one of four modes:
identity, functional bytecode expressions, monotone I-splines, or
piecewise linear interpolation. These transformations are specified in R
and passed to Stan as data so that the likelihood remains deterministic
and reproducible.

### 3.1 Marker-only random effects (no inner id term)

In some studies, we may want **marker-level random effects only**,
without subject-specific marker deviations. This is done by omitting the
inner `( ... | id )` term inside the marker block. In this
configuration, marker-by-id random effects are disabled (internal
standata `Q_idm = 0`), and covariance-style associations that require
those effects (`corr`, `vcov`) are not permitted. Note: assoc =
c(“corr”) or assoc = c(“vcov”) is not allowed when Q_idm = 0.

Code

``` r

# Marker-only random effects: no inner ( ... | id )
formulaLong_marker_only <- y ~ time + x1 +
  (1 + time || id) +
  (0 + x1 || marker)
```

Scales of the parameters can be set via
[`joinme_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md),
which allows explicit control over the prior distributions for fixed
effects, association coefficients, and correlation matrices. Stricter
shrinkage can be set via `shrinkage` params.

Code

``` r

# Define explicit priors for transparency (optional but recommended)
# Note: JoiNMe exposes beta (fixed effects), alpha (associations), and lkj
priors <- joinme_priors(
  beta = list(scale = 2.5),       # Fixed effects (Student-t, df = 6)
  alpha = list(scale = 1.0),      # Association coefficients
  lkj = 2.0                       # Correlation matrices (regularizing)
)
```

## 4 Model Fitting

We fit the model using Hamiltonian Monte Carlo (HMC) via Stan. By
default, `joinme` prefers `cmdstanr` if compiled, otherwise it falls
back to `rstan`. The control parameters can be adjusted using `cmdstanr`
syntax.

Code

``` r

has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL))
engine <- if (has_cmdstan) "cmdstanr" else "rstan"

control <- list(
  engine = engine,
  chains = 1,           # Use 4 chains for publication
  iter_warmup = 200,    # Typically 1000+
  iter_sampling = 200,  # Typically 1000+
  parallel_chains = 1,
  adapt_delta = 0.90,   # Increase if divergences occur
  max_treedepth = 12
)
if (engine == "cmdstanr") control$init <- 0

fit <- joinme(
  dataLong = dataLong,
  dataEvent = dataEvent,
  # Longitudinal formula: Fixed effects + Random effects
  formulaLong = y ~ time + x1 + (1 + time || id) + (0 + x1 + (1 + time || id) || marker),
  # Survival formula: Baseline covariates
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  # Association structure: Link hazard to the total current value of markers
  assoc = c("cv_total"),
  transforms = joinme_tf(cv_total = "identity"),
  priors = priors,
  control = control
)
#> 
#> SAMPLING FOR MODEL 'joinme_fit_threading' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.005868 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 58.68 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 400 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 400 [ 25%]  (Warmup)
#> Chain 1: Iteration: 200 / 400 [ 50%]  (Warmup)
#> Chain 1: Iteration: 201 / 400 [ 50%]  (Sampling)
#> Chain 1: Iteration: 300 / 400 [ 75%]  (Sampling)
#> Chain 1: Iteration: 400 / 400 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 291.331 seconds (Warm-up)
#> Chain 1:                544.145 seconds (Sampling)
#> Chain 1:                835.476 seconds (Total)
#> Chain 1:

# Counting-process fit with left truncation/time-split covariates:
fit_split <- joinme(
  dataLong = sim_split$dataLong,
  dataEvent = dataEvent_split,
  formulaLong = y ~ time + x1 + (1 + time || id) + (0 + x1 + (1 + time || id) || marker),
  formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
  assoc = c("cv_total"),
  control = control
)
#> 
#> SAMPLING FOR MODEL 'joinme_fit_threading' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.008307 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 83.07 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 400 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 400 [ 25%]  (Warmup)
#> Chain 1: Iteration: 200 / 400 [ 50%]  (Warmup)
#> Chain 1: Iteration: 201 / 400 [ 50%]  (Sampling)
#> Chain 1: Iteration: 300 / 400 [ 75%]  (Sampling)
#> Chain 1: Iteration: 400 / 400 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 308.23 seconds (Warm-up)
#> Chain 1:                184.166 seconds (Sampling)
#> Chain 1:                492.396 seconds (Total)
#> Chain 1:
```

## 5 Convergence Diagnostics

We use the `posterior` and `bayesplot` packages to inspect the MCMC
chains.

### 5.1 R-hat and Effective Sample Size

Values of $`\hat{R} < 1.01`$ and ESS \> 400 indicate reasonable
convergence for inference, as recommended by Stan developers.

Code

``` r

# robust summary of regression coefficients on the renamed user-facing scale
draws_obj <- draws(fit, format = "draws_array")

draws_sub <- posterior::subset_draws(
  draws_obj,
  variable = c("^\\(Intercept\\)$", "^time$", "^event", "^cv_total"),
  regex = TRUE
)

posterior::summarise_draws(
  draws_sub,
  default_summary_measures(),
  default_convergence_measures()
)
#> # A tibble: 3 × 10
#>   variable     mean median    sd   mad       q5   q95  rhat ess_bulk ess_tail
#>   <chr>       <dbl>  <dbl> <dbl> <dbl>    <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept) 0.837  0.851 0.426 0.265  0.250   1.33  0.995     264.     184.
#> 2 time        0.177  0.179 0.117 0.108 -0.00460 0.355 1.00      393.     190.
#> 3 cv_total    0.429  0.345 0.389 0.226  0.0712  1.20  0.996     149.     178.
```

### 5.2 Traceplots

Visual inspection of traceplots helps identify mixing issues or stuck
chains. Note that the marker-specific association term for the first
association (if marker weights are shared) is forced to be positive. If
`shared_marker_weights == FALSE`, then the marker-specific association
terms are all constrained to be positive.

Code

``` r

# Plot traces for the association parameter using the helper API
mcmc_plot(fit, variable = "cv_total", type = "trace")
```

![](joinme-workflow_files/figure-html/trace-1.png)

## 6 Dynamic Prediction

One of the most powerful features of joint models is **dynamic
prediction**: updating survival probabilities as new biomarker data
becomes available.

If competing risks are present, the dynamic survival and cumulative
hazard outputs correspond to the overall survival that integrates the
sum of cause-specific hazards. Cause-specific survival curves are not
returned by the default prediction outputs.

We select a subject who had an event at a later time and predict their
survival probability conditionally on data observed up to an earlier
time point.

Code

``` r

# Select a subject with an event > 5
target_id <- dataEvent$id[which(dataEvent$time > 5 & dataEvent$event == 1)[1]]
if (is.na(target_id)) target_id <- dataEvent$id[1]

# Define the conditioning time (time_start)
t_cond <- 4

# You can also pass a column name in newdataEvent, e.g. time_start = "time_start".

# Filter history up to t_cond
ndLong <- dataLong |> filter(id == target_id, time <= t_cond)
ndEvent <- dataEvent |> filter(id == target_id)

print(paste("Predicting for ID:", target_id, "conditioned on history up to t =", t_cond))
#> [1] "Predicting for ID: 4 conditioned on history up to t = 4"
```

Run the prediction:

Code

``` r

# Predict conditional survival from t_cond onwards
preds <- tryCatch(
  predict(
    fit,
    newdataLong = ndLong,
    newdataEvent = ndEvent,
    time_start = t_cond,
    times = seq(t_cond, 8, length.out = 30),
    process = "event", # Focus on survival and cumhaz
    control = list(progress = FALSE) # just to silence the progress bar
  ),
  error = function(e) {
    message("Prediction skipped: ", conditionMessage(e))
    NULL
  }
)

head(preds$predictions$survival)
#>     id     time  Survival    Median  Est.Error       L95       U95
#> 4.1  4 4.000000 1.0000000 1.0000000 0.00000000 1.0000000 1.0000000
#> 4.2  4 4.137931 0.9674819 0.9681881 0.01268423 0.9396159 0.9871276
#> 4.3  4 4.275862 0.9360509 0.9373449 0.02482678 0.8804698 0.9744973
#> 4.4  4 4.413793 0.9056763 0.9075484 0.03643158 0.8226668 0.9621030
#> 4.5  4 4.551724 0.8763169 0.8787884 0.04750484 0.7663613 0.9499354
#> 4.6  4 4.689655 0.8479239 0.8509846 0.05805437 0.7116180 0.9379819
```

`concordance(fit, ...)` now evaluates dynamic prediction on a dense
per-subject survival grid between the landmark and the requested
horizon. This keeps the discrimination calculation aligned with the same
prediction machinery used by the fitted plotting helpers and avoids the
instability that can arise from endpoint-only event prediction.

Finally, we visualise the predicted survival curve. The shaded area
represents the 95% credible interval, accounting for uncertainty in both
the parameters and the random effects.

Code

``` r

if (!is.null(preds)) {
  survival_plot(preds) +
    ggtitle(paste("Dynamic Survival Prediction for Subject", target_id)) +
    geom_vline(xintercept = dataEvent$time[dataEvent$id == target_id], linetype = "dashed", color = "red") +
    annotate("text", x = dataEvent$time[dataEvent$id == target_id], y = 0.1, label = "Observed event time", color = "red", hjust = 1.1)
}
```

![](joinme-workflow_files/figure-html/predict-plot-1.png)

## 7 Posterior reporting helpers

The scientific reporting layer should use the same user-facing parameter
names as the tables returned by `summary(fit)`. The
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md) and
[`mcmc_plot()`](https://rdrr.io/pkg/joinme/man/mcmc_plot.html) helpers
make that possible without directly handling raw Stan variables.

Code

``` r

renamed_draws <- draws(fit, variables = c("time", "cv_total"), format = "draws_df")
head(renamed_draws)
#> # A draws_df: 6 iterations, 1 chains, and 2 variables
#>    time cv_total
#> 1  0.14     0.34
#> 2  0.23     0.61
#> 3  0.13     0.37
#> 4  0.28     0.17
#> 5  0.36     0.21
#> 6 -0.18     0.50
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

mcmc_plot(fit, variable = c("time", "cv_total"), type = "intervals")
```

![](joinme-workflow_files/figure-html/posterior-reporting-1.png)

The plotting wrappers expose the same fitted and predicted views as the
main [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
while making the intent explicit in scripts and reports:

Code

``` r

longitudinal_plot(fit, subject = 1, scale = "epred")
association_plot(fit)
diagnostic_plot(fit, type = "ess_bulk")

if (!is.null(preds)) {
  longitudinal_plot(preds)
  survival_plot(preds)
  cumhaz_plot(preds)
}
```

For fitted longitudinal displays, the heatmap is now an alternative to
the usual separated marker curves:

Code

``` r

plot(
  fit,
  type = "longitudinal",
  longitudinal_style = "heatmap",
  scale = "epred",
  threshold = 0.05
)
```

`plot(fit, ...)` uses fitted posterior samples directly and does not
accept a conditioning argument. Named covariate profiles belong to the
conditional- effects estimand. Build those profiles with
`make_conditions()` and pass them to
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html):

Code

``` r

profiles <- make_conditions(dataEvent, vars = c("x1", "x2"))

effects_data <- conditional_effects(
  fit,
  effects = list(
    longitudinal = "time",
    event = c("x1", "x2")
  ),
  conditions = profiles,
  process = c("longitudinal", "event"),
  longitudinal_estimand = "population",
  plot = FALSE
)

plot(effects_data, ask = FALSE)
```

The longitudinal calculation has three estimands. `"population"` uses
only the fixed-effect contribution. `"marker"` adds the fitted
marker-level deviation while still excluding subject and
marker-by-subject deviations. `"marginal_marker"` calculates each
selected marker response first and averages those responses within every
posterior draw. The last ordering is important for nonlinear links: it
estimates the average expected response rather than the expected
response at an average marker effect. The `markers` argument restricts
the marker-specific panels or the equally weighted set entering that
average.

Code

``` r

average_marker_effect <- conditional_effects(
  fit,
  effects = "time",
  process = "longitudinal",
  longitudinal_estimand = "marginal_marker",
  markers = c("m1", "m2"),
  plot = FALSE
)
```

The event tables contain the covariate-specific relative hazard
`exp(W * gamma)`, holding the baseline hazard and longitudinal
association contribution fixed.
