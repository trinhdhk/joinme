# Prepare standata for Joint Nested Mixed-effects (JoiNMe) model.

This function process the input from `fit` to prepare standata for Joint
Nested Mixed-effects (JoiNMe) model.

## Usage

``` r
joinme_standata(
  formulaLong,
  dataLong,
  formulaEvent,
  dataEvent,
  formulaVCov = ~1,
  formulaDist = NULL,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  eps_fd = 0.01,
  assoc = c("cv_mean"),
  families = NULL,
  transforms = NULL,
  prior_specification = NULL,
  allow_marker_crosscorr = 1L,
  shrinkage = 0L,
  basehaz = c("bs", "ns", "formula"),
  basehaz_n_knots = 5L,
  basehaz_knots = NULL,
  basehaz_degree = 3L,
  basehaz_formula = ~1 + time,
  tau_spline = 0.4,
  quadrature_nodes = NULL,
  vcov_diag_link = c("softplus", "exp"),
  mixture = NULL,
  include_survival = TRUE,
  seed = .Random.seed[[1]]
)
```

## Arguments

- formulaLong:

  Longitudinal formula defining fixed effects, id-level effects, and the
  marker block. The marker block may optionally include an inner
  `( ... | id )` term for marker-by-id random effects. When the inner
  term is omitted, marker-by-id random effects are disabled (Q_idm = 0).
  Use `||` to enforce independence. In nested marker terms, outer
  `( ... || marker )` keeps marker-only and marker-by-id blocks
  independent, while inner `( ... || id )` makes the marker-by-id
  covariance diagonal. Grouping terms may use
  `weighted(group, weights = <column>)` to define formula-scoped
  positive subject/group weights.

- dataLong:

  Long-format longitudinal data with columns for id, marker, time,
  outcome, and covariates referenced in `formulaLong`.

- formulaEvent:

  Survival formula for baseline covariates and event model. Supported
  LHS forms are:

  - `survival::Surv(time, status)`

  - `survival::Surv(start, stop, status)`

  - `survival::Surv(time, status, type = "left")`

  - `survival::Surv(time1, time2, type = "interval2")` The former
    `type = "interval"` form is intentionally rejected. For an
    interval-censored event in \\(L,R\]\\, the prepared data contain an
    event-free row from zero to \\L\\ and a failure-within-interval row
    from \\L\\ to \\R\\. Their log-likelihood contributions sum to
    \\\log\\S(L)-S(R)\\\\. This is a full event likelihood, not a Cox
    partial likelihood.

- dataEvent:

  Event-process data with either one row per id (`Surv(time, status)`)
  or multiple time-split rows per id (`Surv(start, stop, status)`).
  Covariates referenced in `formulaEvent` may vary by interval in the
  split form. Left- and interval-censored responses require one row per
  id and presently describe one event type.

- formulaVCov:

  Covariance regression specification for id-specific marker-by-id
  effects. A formula applies the same observed covariates to both
  components. A named list, `list(sd = ~ ..., corr = ~ ...)`, regresses
  standard deviations and off-diagonal partial correlations
  independently. If the marker block omits the inner `( ... | id )`,
  then marker-by-id effects are absent and covariance-style associations
  (`corr`, `vcov`) are not allowed. Downstream, `corr` uses the
  off-diagonal entries of the subject-specific Cholesky-correlation
  factor `K`, whereas `vcov` uses those off-diagonal `K` entries
  together with the subject-specific standard deviations. If both `corr`
  and `vcov` are requested, `vcov` is kept and `corr` is ignored with a
  warning. The default `~ 1` is valid and yields an intercept-only
  covariance regression with no subject-level slope columns.

- formulaDist:

  Optional list of formulas for distributional regression. Supported
  parameters are `sigma`, `nu`, `phi`, `alpha`, `kappa`, and `tau`. Here
  `nu` is reserved for Student-t degrees of freedom; `kappa` is the
  positive Beta sample-size parameter defining shapes \\\mu\kappa\\ and
  \\(1-\mu)\kappa\\; and `tau` is the skew-double-exponential
  quantile/asymmetry parameter in \\(0,1)\\.

  Three forms are supported:

  1.  Named list with RHS-only formulas, e.g.
      `list(sigma = ~ 1 + time)`.

  2.  Unnamed list with LHS parameter names, e.g.
      `list(sigma ~ 1 + time)`. Random-effects terms with `|` are
      supported; nested random-effects formulas are not. Grouping
      factors in distributional random-effects terms may also use
      `weighted(group, weights = <column>)`.

- id_var:

  Column name for subject id in both longitudinal and event data.

- marker_var:

  Column name for marker/biomarker id in longitudinal data.

- time_var:

  Column name for longitudinal time in `dataLong`.

- eps_fd:

  Positive finite-difference step for association derivatives.

- assoc:

  Character vector specifying association components (e.g., "cv_mean",
  "cs_total", "corr", "vcov"). The `corr` channel targets off-diagonal
  entries of the subject-specific Cholesky-correlation factor `K`; the
  `vcov` channel targets those off-diagonal `K` entries together with
  the subject-specific standard deviations.

- families:

  Optional marker-specific family specification. Can be character family
  names or `jm_family(...)` entries with per-marker links. If NULL, all
  markers use Gaussian responses with identity link. Supported named
  forward links are `identity`, `log`, `logit`, `probit`, and `exp`;
  [`jm_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md)
  also accepts an invertible formula link or a directly specified
  inverse-link formula. In formula syntax, `inv_Phi`/`qnorm`/`probit`
  denote the standard normal quantile and `Phi`/`pnorm` denote the
  standard normal CDF. Thus a probit forward link is inverted to `Phi`
  before Stan data are constructed. A skew-Laplace marker may
  additionally fix its quantile/asymmetry parameter with
  `jm_family("skew_laplace", tau = 0.8)`. Fixed values must lie strictly
  between zero and one and are carried separately for each marker.

- transforms:

  Optional list specifying transformations for association terms. Each
  element (cv_total, cs_total, corr, vcov) is a list with a `type` and
  fields required by that type (see
  [`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md)).
  For covariance-style terms, transforms apply either to off-diagonal
  `K` features (`corr`) or to the combined off-diagonal `K` plus
  subject-specific SD features (`vcov`). Functional transforms may
  request fit-only affine-shift parameters through `intercept = TRUE`
  and/or `slope = TRUE` inside the transform formula. Ordered
  piecewise-linear fits use
  `list(type = "pwlin", knots = ..., direction = "increasing")` (or
  `"decreasing"`). Their knot ordinates are estimated from simplex
  increments; the earlier `x` field is accepted as an alias for `knots`,
  while the earlier `y` field no longer fixes the fitted curve.

- prior_specification:

  A component-based
  [`jm_prior()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md)
  declaration. Global `intercept` and `slope` declarations are inherited
  unless a scientific component supplies the corresponding role
  explicitly. Marker-weighted associations are governed entirely by its
  `marker_weights` component. `marker_weights$offset` is a wholly named
  or wholly unnamed numeric vector; with term-specific weights it may
  instead be a named collection indexed by `cv_total`, `cs_total`,
  `cv_marker`, or `cs_marker`. `marker_weights$family = "constant"` (or
  `"none"`) uses the offset exactly and fits no marker-weight
  coefficient. A stochastic family fits
  `offset + common location + unit-scale marker departure`, with the
  common-location fitting prior in `marker_weights$intercept`. It does
  not supply a population value to simulation. The association slope
  supplies the multiplier for the completed weighted marker feature.

- allow_marker_crosscorr:

  Integer flag; 1 allows cross-marker correlation in marker RE.

- shrinkage:

  Compatibility flag used by latent-class component distributions and
  simulation: 0 = Student-t(6), 1 = Laplace, and 2 = Normal. Ordinary
  coefficient and marker priors are declared separately through
  [`jm_prior()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).

- basehaz:

  Baseline hazard basis type: "bs", "ns", or "formula".

- basehaz_n_knots:

  Number of internal knots for spline baseline hazards.

- basehaz_knots:

  Optional numeric vector of internal knots for spline baseline hazards.

- basehaz_degree:

  Degree of spline basis for baseline hazard.

- basehaz_formula:

  Formula for baseline hazard when `basehaz = "formula"`.

- tau_spline:

  Prior scale for spline coefficients (penalised spline).

- quadrature_nodes:

  Optional positive integer target for total quadrature points. Allowed
  values are exactly `7`, `15`, `31`, `41`, `51`, and `61`.

- vcov_diag_link:

  Link for the subject-specific standard deviation regression:
  "softplus" or "exp".

- mixture:

  Internal latent-progress mixture specification. The public interface
  is
  [`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md);
  ordinary
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md)
  calls leave this as `NULL`, which contributes no mixture parameters or
  mixture prior.

- include_survival:

  Logical indicator that the supplied event rows represent an observed
  survival process. The fitting entry points set this to `FALSE` only
  when they have constructed an internal likelihood-neutral scaffold for
  a longitudinal-only analysis.

- seed:

  Optional random seed for deterministic components of standata.

- flag_resid_dim:

  Integer flag to include residual dimension checks.

## Details

This builder performs three key steps:

1.  Constructs fixed-effect and random-effect design matrices.

2.  Creates distributional regression matrices (including optional
    random effects).

3.  Assembles spline bases and Gauss-Kronrod nodes for survival
    integration.
