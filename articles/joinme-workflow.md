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

1.  **Data assembly (R)**: build design matrices, scale time to \[0,
    1\], and prepare Gauss-Kronrod nodes (default 15, configurable) for
    survival integration.
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
    [`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html),
    and explicit plotting wrappers such as
    [`longitudinal_plot()`](https://trinhdhk.github.io/joinme/reference/longitudinal_plot.html)
    and
    [`survival_plot()`](https://trinhdhk.github.io/joinme/reference/plot_helpers.html).

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
#
# For a probit marker, use link = "probit" or link = ~ inv_Phi(x).
# A direct inverse-link declaration instead uses inv_link = ~ Phi(x).

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
#>   id marker      time         x1        x2         y
#> 1  1     m1 0.0000000 -0.2942841 0.1433918 3.2878746
#> 2  1     m1 0.7272727 -0.2942841 0.1433918 0.4591744
#> 3  1     m1 1.4545455 -0.2942841 0.1433918 1.5627928
#> 4  1     m1 2.1818182 -0.2942841 0.1433918 0.7059543
#> 5  1     m1 2.9090909 -0.2942841 0.1433918 0.6613724
#> 6  1     m1 3.6363636 -0.2942841 0.1433918 4.6048840
head(dataEvent)
#>   id         x1         x2      time event time_start time_stop
#> 1  1 -0.2942841  0.1433918 4.2850991     1          0 4.2850991
#> 2  2 -0.5631947  0.7502855 0.1089635     1          0 0.1089635
#> 3  3  0.3776016 -0.1379869 8.0000000     0          0 8.0000000
#> 4  4 -0.1430829  0.6646221 0.5845890     1          0 0.5845890
#> 5  5  0.2235321 -1.0430229 8.0000000     0          0 8.0000000
#> 6  6  0.4042901  0.5142336 8.0000000     0          0 8.0000000
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
#> 15 25
prop.table(table(dataEvent$event))
#> 
#>     0     1 
#> 0.375 0.625
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
Student-t, Poisson, negative binomial, Bernoulli, beta, or ordinal).
When a family includes additional parameters (for example \sigma, \nu,
or \phi), those parameters can be modelled through `formulaDist`. This
allows heteroscedasticity or covariate-dependent dispersion while
keeping the mean structure aligned across markers.

#### 3.0.2 Association transformations

The association terms can be transformed using one of four modes:
identity, functional bytecode expressions, monotone I-splines, or
ordered piecewise-linear interpolation. These transformations are
specified in R and passed to Stan as data so that the likelihood remains
deterministic and reproducible.

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
  iter_warmup = 100,    # Typically 1000+
  iter_sampling = 100,  # Typically 1000+
  parallel_chains = 1,
  adapt_delta = 0.8,   # Increase if divergences occur
  max_treedepth = 12
)

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
#> Chain 1: Gradient evaluation took 0.005799 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 57.99 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 40.707 seconds (Warm-up)
#> Chain 1:                65.622 seconds (Sampling)
#> Chain 1:                106.329 seconds (Total)
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
#> Chain 1: Gradient evaluation took 0.008789 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 87.89 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 134.044 seconds (Warm-up)
#> Chain 1:                90.486 seconds (Sampling)
#> Chain 1:                224.53 seconds (Total)
#> Chain 1:
```

## 5 Convergence Diagnostics

We use the `posterior` and `bayesplot` packages to inspect the MCMC
chains.

### 5.1 R-hat and Effective Sample Size

Values of \hat{R} \< 1.01 and ESS \> 400 indicate reasonable convergence
for inference, as recommended by Stan developers.

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
#>   variable     mean median    sd    mad      q5   q95  rhat ess_bulk ess_tail
#>   <chr>       <dbl>  <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept) 1.03  1.16   0.575 0.539  0.110   1.83  1.01     102.      57.6
#> 2 time        0.253 0.264  0.111 0.131  0.0839  0.431 0.998     87.9     71.3
#> 3 cv_total    0.149 0.0804 0.234 0.0876 0.00460 0.425 1.01      51.5     76.6
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
#> [1] "Predicting for ID: 7 conditioned on history up to t = 4"
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
#> 7.1  7 4.000000 1.0000000 1.0000000 0.00000000 1.0000000 1.0000000
#> 7.2  7 4.137931 0.9813482 0.9848616 0.01081059 0.9482546 0.9914852
#> 7.3  7 4.275862 0.9629472 0.9697253 0.02141386 0.8966224 0.9830937
#> 7.4  7 4.413793 0.9447983 0.9546153 0.03178993 0.8453193 0.9748782
#> 7.5  7 4.551724 0.9269025 0.9395565 0.04191772 0.7951691 0.9666450
#> 7.6  7 4.689655 0.9092609 0.9246760 0.05177424 0.7461024 0.9583546
```

`concordance(fit)` is a follow-up-wide survival-curve concordance rather
than an AUC at a selected horizon. By default, residual follow-up for
subject (i) starts at their final longitudinal measurement strictly
before the observed event or censoring time, (T\_{0i}). For an event at
residual time (r_i=T_i-T\_{0i}), the event subject is compared with
every subject known to survive longer. The pair is concordant when
(S_i(r_i)\<S_j(r_i)): both conditional survival curves are evaluated at
the same earlier event time. This is the survival-curve concordance
proposed by Antolini and colleagues and does not use the outcome time to
construct a fixed subject score.

Premature censoring is handled through observability. A subject censored
before (r_i) is not comparable with that event; censoring at (r_i)
establishes survival through that time and remains comparable.
`type_weights = "n/G2"` requests Uno’s inverse-censoring weighting
through
[`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html).
Fitted monotone splines and ordered piecewise-linear associations enter
every conditional curve through their posterior ordinates, exactly as in
[`predict()`](https://rdrr.io/r/stats/predict.html).

`auc(fit, ...)` retains the horizon-specific estimand. At landmark (s)
and horizon (t), it compares cumulative risk (1-S_i(ts)) between cases
observed by (t) and controls known to remain event-free beyond (t).
Subjects censored earlier have unknown case/control status and are
omitted. Dynamic prediction for both methods reuses the exact fitted
B-spline, natural-spline, or formula baseline-hazard basis, including
its original centring constants.

The estimands, censoring rules, event-time weighting, computational
steps, and reporting recommendations are derived in the [dedicated
discrimination
vignette](https://trinhdhk.github.io/joinme/articles/joinme-discrimination.md).

Code

``` r

concordance(fit)
concordance(fit, time_start = 2, type_weights = "n/G2")
auc(fit, time_start = 2, time_horizon = 5)
```

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
[`mcmc_plot()`](https://trinhdhk.github.io/joinme/reference/mcmc_plot.html)
helpers make that possible without directly handling raw Stan variables.

Code

``` r

renamed_draws <- draws(fit, variables = c("time", "cv_total"), format = "draws_df")
head(renamed_draws)
#> # A draws_df: 6 iterations, 1 chains, and 2 variables
#>   time cv_total
#> 1 0.32   0.0079
#> 2 0.35   0.0704
#> 3 0.29   0.0444
#> 4 0.37   0.1085
#> 5 0.28   2.0446
#> 6 0.44   0.0792
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
