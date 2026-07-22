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
  beta_prior = NULL,
  alpha_prior = NULL,
  iota_prior = NULL,
  lkj_prior = NULL,
  allow_marker_crosscorr = 1L,
  shrinkage = 0L,
  marker_weights = NULL,
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  basehaz = c("bs", "ns", "formula"),
  basehaz_n_knots = 5L,
  basehaz_knots = NULL,
  basehaz_degree = 3L,
  basehaz_formula = ~1 + time,
  tau_spline = 0.4,
  quadrature_nodes = NULL,
  vcov_diag_link = c("softplus", "exp"),
  tau_fixed = NULL,
  seed = NULL
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

  - `survival::Surv(time1, time2, type = "interval2")` The legacy
    `type = "interval"` form is intentionally rejected.

- dataEvent:

  Event-process data with either one row per id (`Surv(time, status)`)
  or multiple time-split rows per id (`Surv(start, stop, status)`).
  Covariates referenced in `formulaEvent` may vary by interval in the
  split form.

- formulaVCov:

  Covariance regression formula for id-specific marker-by-id effects. If
  the marker block omits the inner `( ... | id )`, then marker-by-id
  effects are absent and covariance-style associations (`corr`, `vcov`)
  are not allowed. Downstream, `corr` uses the off-diagonal entries of
  the subject-specific Cholesky-correlation factor `K`, whereas `vcov`
  uses those off-diagonal `K` entries together with the subject-specific
  standard deviations. If both `corr` and `vcov` are requested, `vcov`
  is kept and `corr` is ignored with a warning. The default `~ 1` is
  valid and yields an intercept-only covariance regression with no
  subject-level slope columns.

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
  markers use Gaussian responses with identity link. Supported
  links/inverse-links: `identity`, `log`, `logit`, `probit`, `exp`.

- transforms:

  Optional list specifying transformations for association terms. Each
  element (cv_total, cs_total, corr, vcov) is a list with a `type` and
  fields required by that type (see
  [`build_standata_transforms()`](https://trinhdhk.github.io/joinme/reference/build_standata_transforms.md)).
  For covariance-style terms, transforms apply either to off-diagonal
  `K` features (`corr`) or to the combined off-diagonal `K` plus
  subject-specific SD features (`vcov`). Functional transforms may
  request fit-only affine-shift parameters through `intercept = TRUE`
  and/or `slope = TRUE` inside the transform formula.

- beta_prior:

  Prior specification for longitudinal fixed effects.

- alpha_prior:

  Prior specification for association parameters.

- iota_prior:

  Prior specification for fit-only affine-shift intercept and slope
  parameters used by functional association transforms.

- lkj_prior:

  Prior specification for correlation structures.

- allow_marker_crosscorr:

  Integer flag; 1 allows cross-marker correlation in marker RE.

- shrinkage:

  Integer flag controlling shrinkage behaviour for marker-by-id effects.
  0 = student_t(6, 0, 1), 1 = double_exponential(0, 1), 2 =
  std_normal();

- marker_weights:

  Optional base weights for marker-specific association components. If
  `shared_marker_weights = TRUE`, provide one numeric vector of length D
  (or one named numeric vector keyed by marker level) that is shared
  across all active weighted association terms. If
  `shared_marker_weights = FALSE`, you may instead provide a named list
  with entries `cv_total`, `cs_total`, `cv_marker`, and `cs_marker`.
  Each supplied entry is aligned to marker order and used only for the
  matching active weighted association term.

- fixed_marker_weights:

  Logical; if TRUE, marker weights are kept fixed at `marker_weights`
  (no perturbation). If FALSE, marker weights are estimated via signed
  additive perturbations, `marker_weights + z_marker_weights`. The
  standardized latent family is selected by `shrinkage`: Student-t(6),
  Laplace, or Normal for 0, 1, or 2.

- shared_marker_weights:

  Logical; if TRUE, all active weighted marker-based association terms
  share one marker-weight structure. If FALSE, each active weighted
  marker-based association term gets its own marker-weight structure.

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
  values are exactly `7`, `15`, `31`, `41`, `51`, and `61`. Only the
  node count is passed to Stan; GK nodes/weights are fixed in Stan.

- vcov_diag_link:

  Link for the subject-specific standard deviation regression:
  "softplus" or "exp".

- tau_fixed:

  Optional fixed quantile/asymmetry parameter for the skew double
  exponential family. It must lie strictly between zero and one; the
  value \\0.5\\ gives the symmetric double exponential distribution. It
  cannot be combined with a `tau` distributional regression.

- seed:

  Optional random seed for deterministic components of standata.

- flag_resid_dim:

  Integer flag to include residual dimension checks.

## Details

This builder performs three key steps:

1.  Constructs fixed-effect and random-effect design matrices on a
    scaled time axis.

2.  Creates distributional regression matrices (including optional
    random effects).

3.  Assembles spline bases and Gauss-Kronrod nodes for survival
    integration. Only the node count is passed to Stan; nodes/weights
    are hardcoded in the Stan functions that evaluate the cumulative
    hazard.

Time is scaled internally as `t_scaled = t/t_max` for numerical
stability; indices are recorded to rescale time-associated coefficients
back to the original units inside Stan.
