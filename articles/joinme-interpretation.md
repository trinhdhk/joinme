# Conditional effects and contrasts in JoiNMe

## Purpose

This vignette defines the fitted quantities returned by
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
and
[`conditional_contrast()`](https://trinhdhk.github.io/joinme/reference/conditional_contrast.md)
and demonstrates their public interface. It also states the
event-likelihood quantities required to distinguish interval censoring,
delayed entry and ordinary right censoring.

The examples use generic variable names that can be adapted to a study.

## The event model is a full likelihood

Let T_i be the event time for subject i. Conditional on the fitted
random effects and model parameters, the cause-specific hazard can be
formulated as

h\_{ik}(t) =h\_{0k}(t) \exp\\\boldsymbol w_i(t)^\top\boldsymbol\gamma_k
+a_i(t)\\,

where h\_{0k}(t) is the baseline hazard for cause k, \boldsymbol
w_i(t)^\top\boldsymbol\gamma_k is the direct event-regression
contribution, and a_i(t) contains the selected longitudinal association
terms. The all-cause cumulative hazard and survival function are

H_i(t)=\int_0^t\sum_k h\_{ik}(u)\\du, \qquad S_i(t)=\exp\\-H_i(t)\\.

`joinme` estimates the baseline hazard and evaluates these integrals by
Gauss–Kronrod quadrature. The event contribution is therefore a full
event-time likelihood. It is not the Cox partial likelihood. This
distinction matters for left and interval censoring: their probabilities
depend on the absolute survival curve, including the baseline hazard,
and cannot be obtained from event ordering alone.

### What each censoring statement contributes

For one event type, the four supported observations have the following
meaning. Endpoint conventions follow
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html).

| Observation           | Knowledge about T_i | Likelihood contribution |
|-----------------------|--------------------:|------------------------:|
| Right censored at L_i |            T_i\>L_i |                S_i(L_i) |
| Exact event at t_i    |             T_i=t_i |        h_i(t_i)S_i(t_i) |
| Left censored at R_i  |         T_i\leq R_i |              1-S_i(R_i) |
| Interval censored     |    L_i\<T_i\leq R_i |       S_i(L_i)-S_i(R_i) |

For a proper interval (L_i,R_i\], the fitted calculation uses the exact
identity

\begin{aligned} S_i(L_i)-S_i(R_i) &=S_i(L_i)\left\[1-
\exp\\-H_i(R_i)+H_i(L_i)\\\right\],\\ \log\\S_i(L_i)-S_i(R_i)\\
&=-H_i(L_i)+ \log\left\[1-\exp\\-\[H_i(R_i)-H_i(L_i)\]\\\right\].
\end{aligned}

The data preparation represents these two terms as adjacent risk
intervals: a right-censored interval from zero to L_i, followed by a
failure-within-interval contribution from L_i to R_i. The first term is
essential. Omitting it would condition on survival to L_i, giving a
different likelihood and treating the lower inspection time as if it
were an entry time.

The same calculation is used in
[`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md) and
[`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md).
The latent-class model changes the distribution of the selected
random-effect block; it does not replace the event likelihood.

### Supplying `interval2` data

Use one event row per subject. The first time is the lower inspection
limit and the second time is the upper inspection limit.

``` r

dataEvent_interval <- data.frame(
  id = c(1, 2, 3, 4),
  lower = c(4, 3, NA, 2),
  upper = c(NA, 3, 5, 6),
  treatment = c("control", "active", "control", "active"),
  age = c(62, 58, 71, 66)
)

survival::Surv(
  dataEvent_interval$lower,
  dataEvent_interval$upper,
  type = "interval2"
)
```

These rows mean:

- subject 1 is right censored after time 4 (`upper = NA`);
- subject 2 has an exact event at time 3 (`lower = upper`);
- subject 3 is left censored at time 5 (`lower = NA`);
- subject 4 has an event in ((2,6\]).

The model is fitted by placing the same response on the left of
`formulaEvent`:

``` r

fit_interval <- joinme(
  formulaLong = outcome ~ time + treatment +
    (1 + time | id) +
    (1 + time + (1 + time | id) | marker),
  dataLong = dataLong,
  formulaEvent = survival::Surv(lower, upper, type = "interval2") ~
    treatment + age,
  dataEvent = dataEvent_interval,
  assoc = c("cv_mean", "cv_marker"),
  priors = jm_prior(
    longitudinal = list(slope = prior_normal(0, 1)),
    survival = list(slope = prior_normal(0, 1)),
    assoc = list(slope = prior_normal(0, 1))
  )
)
```

The event-data boundary has the following properties:

- An `interval2` response is a single censoring statement for each
  subject.Repeated event rows for the same subject are rejected.
  Counting-process data describe another observation scheme.
- Event covariates in `formulaEvent` are treated as the subject’s
  covariate profile throughout the represented risk time. The latent
  longitudinal trajectory and its association with the hazard may still
  vary continuously with time.
- The lower limit of an interval-censored observation is an inspection
  time, not a delayed-entry time. Delayed entry is represented
  separately by `Surv(start, stop, status)`.
- Left and interval censoring have no observed event cause. Their
  present likelihood is therefore for one event type. Exact and
  right-censored datamay use the package’s cause-specific competing-risk
  representation.
- `concordance()`,
  [`tvROC()`](https://trinhdhk.github.io/joinme/reference/tvROC.md) and
  [`tvAUC()`](https://trinhdhk.github.io/joinme/reference/tvROC.md)
  require an event time that is known exactly or known to exceed a
  censoring time. They deliberately reject left-and interval-censored
  outcomes because those measures need a different censoring-specific
  estimator.

The fitted `log_lik_surv[i]` is the complete subject-level event
contribution. For an interval-censored subject it includes both survival
to the lower limit and failure before the upper limit. `cumhaz_event[i]`
is the cumulative hazard through the represented risk intervals, and
`surv_prob_event[i]` is its corresponding survival probability at the
final upper limit. The latter is a model quantity at that limit; it is
not the probability of the observed interval itself.

## Model components used by conditional effects

The fitted model contains distinct longitudinal, event and association
components.

### Longitudinal coefficients

The population coefficients in `fixef(fit)$formulaLong` describe the
longitudinal linear predictor. With an identity link they are changes on
the response scale. With a nonlinear link they are changes on the link
scale; the corresponding expected-response trajectories are available
from the prediction methods.

The subject, marker, and marker-by-subject terms are deviations from
that population trajectory. They describe variation in the fitted
population; they are not additional population slopes.
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) reports
these deviations, whereas [`coef()`](https://rdrr.io/r/stats/coef.html)
adds a deviation to its relevant population contribution.

### Event coefficients

An event coefficient \gamma_j is a log hazard ratio conditional on the
other event covariates, the longitudinal association, the random effects
and the baseline hazard. Its exponential is a hazard ratio. This is an
instantaneous comparison, not a ratio of event probabilities over a
chosen follow-up period.

### Association coefficients

An association coefficient describes the change in log hazard for a
one-unit change in its transformed longitudinal summary, conditional on
the rest of the joint model. The unit is determined by the chosen
association and its transformation. For example, a coefficient attached
to `cv_mean` acts on the subject-level current value after the declared
transformation. A coefficient attached to a marker-weighted term acts on
the completed weighted mean, not on one marker in isolation.

If an association transformation is nonlinear, one coefficient alone
does not summarise its shape.
[`posterior_assoc()`](https://trinhdhk.github.io/joinme/reference/assoc.md)
and
[`association_plot()`](https://trinhdhk.github.io/joinme/reference/association_plot.md)
return the posterior association curve over a chosen covariate range.

## Conditional effects

[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
describes fitted outcomes while selected predictors vary and the
remaining predictors are held at declared values. It is useful for
displaying nonlinear terms, interactions and response-scale
trajectories. It does not itself compare two groups.

``` r

profiles <- data.frame(
  age = c(55, 70),
  sex = factor(c("female", "female"), levels = levels(dataEvent$sex)),
  cond__ = c("Age 55", "Age 70")
)

effects <- conditional_effects(
  fit,
  effects = list(
    longitudinal = "time:treatment",
    event = "treatment"
  ),
  conditions = profiles,
  int_conditions = list(
    treatment = c("control", "active")
  ),
  process = c("longitudinal", "event"),
  method = "posterior_epred",
  longitudinal_estimand = "marginal_marker",
  event_scale = "hazard_ratio",
  prob = 0.95,
  robust = TRUE,
  plot = FALSE
)

plot(effects, ask = FALSE)
```

When several effects or processes are requested, they may be retained as
separate figures or placed in one display. Patchwork is used for joint
displays because it preserves each process as a complete `ggplot` with
its own outcome scale. An automatic grid, one row, one column, or a
declared design may be used:

``` r

plot(effects, arrange = "grid", ncol = 2, guides = "collect")
plot(effects, arrange = "row")
plot(effects, arrange = "column")

# The letters refer to panels in their returned order: longitudinal effects
# first, followed by event effects. A hash marks unused space.
plot(effects, arrange = "design", design = "AA\nBB")
```

`widths` and `heights` give relative panel dimensions. With
`arrange = "separate"`, which is the default, `plot(..., plot = FALSE)`
returns the named process-and-effect list. With any joint arrangement it
returns one patchwork object, allowing a common title or another
patchwork annotation to be added before printing.

Each row of `conditions` is a common scientific profile. A numeric
predictor not supplied there is held at its observed mean; a categorical
predictor is held at its first fitted level. Supplying all important
adjustment variables is preferable because it makes the comparison
reproducible and avoids relying on an arbitrary reference level.

`effects` may be a character vector shared across processes, or a named
list with `longitudinal` and `event` entries. `int_conditions` supplies
the values of a predictor that is being varied. One- and two-predictor
effects are supported; two continuous predictors may be shown as a
surface.

### Choosing the longitudinal estimand

The three longitudinal estimands answer different questions.

| Value | Random-effect contribution | Result |
|----|----|----|
| `"population"` | Population coefficients only | One trajectory per condition, averaged over selected marker response functions |
| `"marker"` | Population plus fitted marker deviation | A separate trajectory for every selected marker |
| `"marginal_marker"` | Population plus fitted marker deviation | One trajectory per condition, averaged over selected markers within each posterior draw |

For `method = "posterior_epred"`, each marker’s inverse link is applied
before any marker average. This matters when markers have different
families or nonlinear links. The result is an average of expected marker
responses, not the inverse link of an average linear predictor.
`method = "posterior_linpred"` keeps the comparison on the longitudinal
linear-predictor scale.

The `markers` argument restricts the markers shown by `"marker"` or
averaged by the other two estimands:

``` r

average_response <- conditional_effects(
  fit,
  effects = "time",
  conditions = profiles,
  process = "longitudinal",
  method = "posterior_epred",
  longitudinal_estimand = "marginal_marker",
  markers = c("region_1", "region_2"),
  plot = FALSE
)
```

### Population and fitted-subject estimands

By default, subject and marker-by-subject deviations are excluded. The
result therefore describes the fitted population or marker distribution
rather than one observed subject.

To describe a subject represented during fitting, include the fitted
identifier in `conditions` and request its fitted random effects:

``` r

subject_effect <- conditional_effects(
  fit,
  effects = "time",
  conditions = data.frame(
    id = 17,
    treatment = "active",
    age = 62,
    cond__ = "Subject 17"
  ),
  process = "longitudinal",
  longitudinal_estimand = "marker",
  reuse_fitted_re = TRUE,
  plot = FALSE
)
```

This is a fitted-subject conditional trajectory. Even when
`longitudinal_estimand = "population"`, asking to reuse fitted random
effects adds the selected subject and marker-by-subject deviations, so
the result is no longer a purely population trajectory. An identifier
absent from the fitting data is rejected. New-subject inference requires
dynamic prediction rather than borrowing another subject’s fitted
effects.

### Event conditional effects

The event result reports \exp\\\boldsymbol w^\top\boldsymbol\gamma_k\\
for `event_scale = "hazard_ratio"`, or \boldsymbol
w^\top\boldsymbol\gamma_k for `"log_hazard_ratio"`. It isolates the
direct covariate part of `formulaEvent`. The baseline hazard and
longitudinal association are held fixed. Consequently this result is not
a survival probability and is not a dynamic prediction. Use
`predict(..., process = "event")` when the scientific question concerns
survival beyond a landmark time given observed marker history.

With competing risks, event conditional effects are returned separately
foreach event type. Interval-censored fits presently have one event type
because the cause is not observed within the interval.

### Posterior summaries and retained draws

With `summary = TRUE`, each table contains a posterior centre, posterior
standard error and equal-tailed credible limits. `robust = TRUE` uses
the posterior median as the centre; `robust = FALSE` uses the posterior
mean.

With `summary = FALSE`, the result retains one row per posterior draw
and estimand. `.draw` identifies the paired posterior draw, `.value` is
the fitted quantity, and the remaining columns describe the condition,
marker and event type. This form is suitable for probabilities such as
P(\Delta\>0\mid\text{data}) or for another scientifically defined
summary:

``` r

effect_draws <- conditional_effects(
  fit,
  effects = "time",
  conditions = profiles,
  process = "longitudinal",
  longitudinal_estimand = "marginal_marker",
  summary = FALSE,
  plot = FALSE
)

draw_table <- effect_draws$longitudinal$time
```

## Conditional contrasts

[`conditional_contrast()`](https://trinhdhk.github.io/joinme/reference/conditional_contrast.md)
compares group A with group B within every posterior draw. It then
summarises those paired contrasts. This preserves posterior dependence
between the two predictions; subtracting two printed conditional effects
does not.

`groupA` and `groupB` are named lists or named vectors. Lists are safest
when the profile combines numbers, logical values and factor levels.
Values in the groups override the same variables in `conditions`. All
other values are shared, so the two profiles differ only where declared.

``` r

follow_up_profiles <- data.frame(
  time = seq(0, 24, length.out = 61),
  age = 62,
  sex = factor("female", levels = levels(dataEvent$sex)),
  cond__ = paste0("Month ", seq(0, 24, length.out = 61))
)

treatment_contrast <- conditional_contrast(
  fit,
  groupA = list(treatment = "active"),
  groupB = list(treatment = "control"),
  conditions = follow_up_profiles,
  process = c("longitudinal", "event"),
  method = "posterior_epred",
  longitudinal_estimand = "marginal_marker",
  event_scale = "hazard_ratio",
  prob = 0.95,
  robust = TRUE,
  plot = FALSE
)

plot(
  treatment_contrast,
  condition_variable = "time",
  arrange = "row",
  guides = "collect"
)
```

The same `arrange`, `ncol`, `nrow`, `design`, `widths`, `heights`, and
`guides` arguments are shared by conditional-effect and
conditional-contrast plots. A row is often clear for two processes; a
column is useful when the horizontal axis should occupy more space.

The longitudinal contrast is

\Delta_i(t)=m_A(t)-m_B(t)

on the scale selected by `method`. A positive expected-response contrast
means the fitted expected longitudinal outcome is higher under profile
A. Whether a higher outcome is favourable depends on the scientific
meaning of the marker.

The event contrast is

\operatorname{HR}\_{A:B,k} =\exp\\\eta\_{A,k}-\eta\_{B,k}\\

or its logarithm. A hazard ratio below one means that the direct event
covariate predictor gives a lower instantaneous hazard for A than B,
holding the baseline hazard and longitudinal association common. It does
not by itself give an absolute risk difference.

### A fitted-subject contrast

The same subject effects may be included on both sides of a contrast.
This defines a fitted prediction for one subject under two covariate
profiles.

``` r

subject_contrast <- conditional_contrast(
  fit,
  groupA = list(treatment = "active"),
  groupB = list(treatment = "control"),
  conditions = data.frame(
    id = 17,
    time = seq(0, 24, length.out = 61),
    age = 62,
    sex = factor("female", levels = levels(dataEvent$sex))
  ),
  process = "longitudinal",
  longitudinal_estimand = "marker",
  reuse_fitted_re = TRUE,
  plot = FALSE
)
```

The identifier cannot differ between group A and group B. The
calculation is conditional on the same fitted random effects on both
sides and does not fit a new subject model.

### Posterior probability of a contrast

Draw-level contrasts make posterior probability statements transparent:

``` r

contrast_draws <- conditional_contrast(
  fit,
  groupA = list(treatment = "active"),
  groupB = list(treatment = "control"),
  conditions = follow_up_profiles,
  process = "longitudinal",
  longitudinal_estimand = "marginal_marker",
  summary = FALSE,
  plot = FALSE
)$longitudinal

probability_positive <- aggregate(
  contrast_draws$.value > 0,
  by = list(time = contrast_draws$time),
  FUN = mean
)
names(probability_positive)[2] <- "posterior_probability_above_zero"
```
