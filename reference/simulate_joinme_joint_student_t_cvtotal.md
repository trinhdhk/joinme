# Simulate joint data (Student-t longitudinal; CV_total survival association)

Simulate joint data (Student-t longitudinal; CV_total survival
association)

## Usage

``` r
simulate_joinme_joint_student_t_cvtotal(
  n_id = 50,
  D = 10,
  n_t = 10,
  seed = .Random.seed[[1]],
  include_marker_only = TRUE,
  Q_idm = 2,
  R_id = 2,
  R_mk = 1,
  beta = c(`(Intercept)` = 1, time = 0.5, x1 = 0.4),
  gamma_w = c(x1 = 0.3, x2 = -0.2),
  alpha_cv_total = 0.6,
  h0 = weibull_h0(shape = 1.4, scale = 6),
  time_cens = 8,
  nu = 5,
  sigma_y = 0.4,
  sigma_fun = NULL,
  formulaDist = NULL,
  beta_sigma = NULL,
  beta_nu = NULL,
  integration_control = list(rel.tol = 1e-06, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1, t_max = 50, expand = 1.7, max_expand = 60L)
)
```

## Arguments

- n_id:

  Number of subjects.

- D:

  Number of markers.

- n_t:

  Measurements per (id, marker).

- seed:

  RNG seed.

- include_marker_only:

  Include marker-only random effects block.

- Q_idm:

  Dimension of marker-by-id basis.

- R_id:

  Dimension of id-level basis.

- R_mk:

  Dimension of marker-only basis.

- beta:

  Fixed effects (Intercept, time, x1) on original scale.

- gamma_w:

  Hazard covariates (e.g. x1, x2).

- alpha_cv_total:

  Association coefficient for CV_total.

- h0:

  Baseline hazard function h0(t).

- time_cens:

  Administrative censoring time.

- nu:

  Student-t df.

- sigma_y:

  Student-t scale.

- sigma_fun:

  Optional function to generate observation-level scales. Signature:
  function(id, times, dataEvent) -\> numeric vector.

- formulaDist:

  Optional distributional regression formulas.

- beta_sigma:

  Optional coefficients for sigma regression.

- beta_nu:

  Optional coefficients for nu regression.

- integration_control:

  list passed to integrate().

- root_control:

  list controlling bracketing.

## Value

list(dataLong, dataEvent, truth, helpers, tmax)
