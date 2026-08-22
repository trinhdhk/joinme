# Declare population truths for JoiNMe simulation

`jm_truth()` describes the population quantities used once to generate a
complete data set. Its coefficient hierarchy follows
[`jm_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md):
global `intercept` and `slope` declarations provide fallbacks, whilst
named model components may replace either role.

A finite numeric vector fixes the corresponding population coefficient.
A `prior_*()` declaration draws one coefficient vector once at the
beginning of a simulation. The realised vector is then held constant for
every subject, marker and observation and is recorded in the returned
simulation's `truth`. Thus, a distribution supplied here describes
variation between simulated data sets, not variation between
observations within one data set.

The declaration also contains the baseline hazard and the covariance
parameters for ordinary random-effect blocks. If a random-effect
standard deviation is omitted, one value per coefficient is drawn from
`Exponential(1)`. If its correlation matrix is omitted, one matrix is
drawn from the LKJ distribution selected by `lkj`. These are the same
population distributions used by the Stan fitting model.

## Usage

``` r
joinme_truth(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc_coef = NULL,
  functional = NULL,
  marker = NULL,
  lkj = prior_lkj(1),
  class = NULL,
  basehaz = list(type = "weibull", shape = 1.4, scale = 6),
  re_params = list(id = list(sd = NULL, corr = NULL), marker = list(sd = NULL, corr =
    NULL), dist = list()),
  .validate = TRUE
)

jm_truth(
  ...,
  intercept = NULL,
  slope = NULL,
  longitudinal = NULL,
  survival = NULL,
  baseline = NULL,
  vcov = NULL,
  marker_weights = NULL,
  assoc_coef = NULL,
  functional = NULL,
  marker = NULL,
  lkj = prior_lkj(1),
  class = NULL,
  basehaz = list(type = "weibull", shape = 1.4, scale = 6),
  re_params = list(id = list(sd = NULL, corr = NULL), marker = list(sd = NULL, corr =
    NULL), dist = list()),
  .validate = TRUE
)
```

## Arguments

- ...:

  Named distributional-parameter truths corresponding to the left-hand
  sides of `formulaDist`, including family or marker selectors.

- intercept, slope:

  Global fixed values or `prior_*()` generating distributions for
  intercept and non-intercept population coefficients.

- longitudinal:

  Fixed values or generating distributions for the longitudinal
  population regression. A named list may contain `intercept` and
  `slope` separately.

- survival:

  Fixed values or a generating distribution for event-model covariate
  slopes.

- baseline:

  Fixed values or generating distributions for coefficients of a formula
  baseline hazard.

- vcov:

  Covariance-regression truth. It may be shared between the
  marginal-standard-deviation and correlation regressions, or supplied
  as `list(sd = ..., corr = ...)`. Each part accepts `intercept`,
  `slope`, and `latent` declarations.

- marker_weights:

  A named list containing `offset`, `intercept`, `family`, and `shared`.
  `offset` is a fixed marker-specific contribution; `intercept` is a
  fixed or once-drawn common marker-weight location; `family` governs
  centred unit-scale marker departures; and `shared` determines whether
  weighted association terms share one set.

- assoc_coef:

  Fixed values or a generating distribution for the active
  longitudinal–event association coefficients, in their fitted order.

- functional:

  Fixed values or generating distributions for affine transformation
  intercepts and slopes.

- marker:

  A `prior_*()` family declaration for standardised marker-level
  effects. Its location and ordinary scale must remain zero and one.

- lkj:

  An object created by
  [`prior_lkj()`](https://trinhdhk.github.io/joinme/reference/prior_lkj.md).
  Its concentration governs each random-effect correlation matrix drawn
  when `re_params` omits `corr`.

- class:

  A named list containing `baseline_prob`, `slope`, and `family`,
  following the class roles accepted by
  [`jm_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
  The slope may be fixed or drawn once for each simulation. The family
  is a probability declaration governing every standardised latent
  coordinate drawn conditionally on its class; it is not itself a
  population coefficient draw.

- basehaz:

  Baseline-hazard truth. Supply a hazard function, a formula in `time`,
  a character family name, or a named list such as
  `list(type = "weibull", shape = 1.4, scale = 6)`.

- re_params:

  Named random-effect covariance declarations for `id`, `marker`, and
  distributional regressions under `dist`. Each ordinary block accepts
  `sd` and `corr`; either may be omitted and drawn once from the fitted
  model's corresponding population distribution.

- .validate:

  Logical; if `TRUE`, check the declaration immediately.

## Value

An object of class `joinme_truth` for the `truth` argument of
[`simulate_joinme()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme.md)
or
[`simulate_joinme_mix()`](https://trinhdhk.github.io/joinme/reference/simulate_joinme_mix.md).

## Examples

``` r
generating_truth <- jm_truth(
  longitudinal = list(
    intercept = 0,
    slope = prior_normal(mu = 0.5, scale = 0.2)
  ),
  survival = c(treatment = -0.4),
  assoc_coef = c(cv_mean = 0.3),
  basehaz = list(type = "weibull", shape = 1.2, scale = 7),
  re_params = list(id = list(sd = NULL, corr = NULL)),
  lkj = prior_lkj(2)
)
```
