#' @name joinme_simulation
#' @title Simulation for joinme
#'
#' @importFrom stats rnorm rexp runif plogis rbeta rbinom rnbinom rpois rt integrate uniroot
#'
#' @description
#' Simulator for the joinme joint model aligned with the v11 Stan semantics.
#'
#' Key features:
#' - Multivariate longitudinal outcomes in long format with irregular times.
#' - Survival hazard depends on CV_total(t) = CV_mean(t) + CV_marker(t).
#' - Survival times drawn by inverse transform sampling using:
#'   * numeric integration (integrate) for cumulative hazard
#'   * numeric root solving (uniroot) for event time.
#'
#' @keywords internal
NULL

# File overview:
# - Simulate longitudinal + survival data consistent with Stan semantics.
# - Provide helper utilities for hazards and root finding.

#' Stable softplus
#' @keywords internal
.softplus <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' Weibull baseline hazard factory
#'
#' @param shape Weibull shape (>0).
#' @param scale Weibull scale (>0).
#' @return function h0(t)
#' @export
weibull_h0 <- function(shape = 1.4, scale = 6.0) {
  force(shape)
  force(scale)
  function(t) {
    t <- pmax(t, 0)
    (shape / scale) * (t / scale)^(shape - 1)
  }
}

#' Find a bracketing interval for root solving
#' @keywords internal
.find_bracket <- function(f, lower = 0, upper = 1, upper_max = 50, expand = 1.7, max_expand = 60L) {
  f_lower <- f(lower)
  if (!is.finite(f_lower) || f_lower > 0) {
    return(NULL)
  }
  f_upper <- f(upper)
  if (!is.finite(f_upper)) {
    return(NULL)
  }
  it <- 0L
  while (f_upper < 0 && upper < upper_max && it < max_expand) {
    upper <- min(upper_max, upper * expand)
    f_upper <- f(upper)
    if (!is.finite(f_upper)) {
      return(NULL)
    }
    it <- it + 1L
  }
  if (f_upper >= 0) list(lower = lower, upper = upper) else NULL
}

#' Simulate joint data (Student-t longitudinal; CV_total survival association)
#'
#' @param n_id Number of subjects.
#' @param D Number of markers.
#' @param n_t Measurements per (id, marker).
#' @param seed RNG seed.
#' @param include_marker_only Include marker-only random effects block.
#' @param Q_idm Dimension of marker-by-id basis.
#' @param R_id Dimension of id-level basis.
#' @param R_mk Dimension of marker-only basis.
#' @param beta Fixed effects (Intercept, time, x1) on original scale.
#' @param gamma_w Hazard covariates (e.g. x1, x2).
#' @param alpha_cv_total Association coefficient for CV_total.
#' @param h0 Baseline hazard function h0(t).
#' @param time_cens Administrative censoring time.
#' @param nu Student-t df.
#' @param sigma_y Student-t scale.
#' @param sigma_fun Optional function to generate observation-level scales.
#'   Signature: function(id, times, dataEvent) -> numeric vector.
#' @param formulaDist Optional distributional regression formulas.
#' @param beta_sigma Optional coefficients for sigma regression.
#' @param beta_nu Optional coefficients for nu regression.
#' @param integration_control list passed to integrate().
#' @param root_control list controlling bracketing.
#'
#' @return list(dataLong, dataEvent, truth, helpers, tmax)
#' @export
simulate_joinme_joint_student_t_cvtotal <- function(
  n_id = 50,
  D = 10,
  n_t = 10,
  seed = 42,
  include_marker_only = TRUE,
  Q_idm = 2,
  R_id = 2,
  R_mk = 1,
  beta = c("(Intercept)" = 1.0, "time" = 0.5, "x1" = 0.4),
  gamma_w = c("x1" = 0.3, "x2" = -0.2),
  alpha_cv_total = 0.6,
  h0 = weibull_h0(shape = 1.4, scale = 6.0),
  time_cens = 8.0,
  nu = 5,
  sigma_y = 0.4,
  sigma_fun = NULL,
  formulaDist = NULL,
  beta_sigma = NULL,
  beta_nu = NULL,
  integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1.0, t_max = 50.0, expand = 1.7, max_expand = 60L)
) {
  # Workflow: simulate covariates -> random effects -> longitudinal -> survival
  # - matches Stan semantics for CV_total association
  set.seed(seed)
  # Step: configure optional parallel execution for simulation.

  dist_formulas <- .normalize_formula_dist(formulaDist)
  .validate_dist_formula_scopes(dist_formulas, "student_t")

  # Event-level data (exogenous)
  dataEvent <- data.frame(
    id = seq_len(n_id),
    x1 = rnorm(n_id),
    x2 = rnorm(n_id),
    stringsAsFactors = FALSE
  )
  dataEvent$Xcov1 <- 0.0

  # ---- id RE u_id
  tau_u <- rep(0.6, R_id)
  if (R_id >= 2) tau_u[2] <- 0.4
  Corr_u <- diag(R_id)
  if (R_id == 2) Corr_u <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  Lcorr_u <- t(chol(Corr_u))
  z_u <- matrix(rnorm(n_id * R_id), n_id, R_id)
  L_u <- diag(tau_u, R_id, R_id) %*% Lcorr_u
  u_id <- z_u %*% t(L_u)

  # ---- marker-only v_marker
  if (include_marker_only && R_mk > 0) {
    tau_v <- rep(0.4, R_mk)
    Corr_v <- diag(R_mk)
    Lcorr_v <- t(chol(Corr_v))
    z_v <- matrix(rnorm(D * R_mk), D, R_mk)
    L_v <- diag(tau_v, R_mk, R_mk) %*% Lcorr_v
    v_marker <- z_v %*% t(L_v)
  } else {
    tau_v <- numeric(0)
    Lcorr_v <- matrix(0, 0, 0)
    z_v <- matrix(0, D, 0)
    v_marker <- matrix(0, D, 0)
  }

  # ---- marker-by-id latent z_w (id-specific, iid standard normal)
  z_w_lat <- array(rnorm(n_id * D * Q_idm), dim = c(n_id, D, Q_idm))
  z_w <- z_w_lat

  # ---- id-specific covariance L_i (diagonal)
  alpha_L <- rep(-0.2, Q_idm)
  lambda_L <- rep(0.3, Q_idm)
  z_L <- matrix(rnorm(n_id * Q_idm), nrow = n_id, ncol = Q_idm)

  L_i <- array(0.0, dim = c(n_id, Q_idm, Q_idm))
  for (i in seq_len(n_id)) {
    for (q in seq_len(Q_idm)) {
      lp <- alpha_L[q] + lambda_L[q] * z_L[i, q]
      L_i[i, q, q] <- .softplus(lp)
    }
  }

  # w_idscaled[i,d] = L_i[i] %*% z_w[i,d]
  w_idscaled <- array(0.0, dim = c(n_id, D, Q_idm))
  for (i in seq_len(n_id)) {
    Li <- matrix(L_i[i, , ], Q_idm, Q_idm)
    for (d in seq_len(D)) w_idscaled[i, d, ] <- as.numeric(Li %*% z_w[i, d, ])
  }

  # ---- design rows (original time)
  fixed_row <- function(i, t) c("(Intercept)" = 1, time = t, x1 = dataEvent$x1[i])
  id_row <- function(t) if (R_id == 1) c(1) else c(1, t)
  mk_row <- function(i, t) if (!include_marker_only || R_mk <= 0) numeric(0) else c(dataEvent$x1[i])
  idm_row <- function(t) if (Q_idm == 1) c(1) else c(1, t)
  haz_row <- function(i) c(x1 = dataEvent$x1[i], x2 = dataEvent$x2[i])

  # CV components (vectorized over t)
  cv_mean_i <- function(i, t) {
    # fixed_part: (1, t, x1) * beta
    # id_part: (1, t) * u_id[i,]
    x1_val <- dataEvent$x1[i]
    fixed_part <- beta[1] + beta[2] * t + beta[3] * x1_val
    id_part <- if (R_id == 1) u_id[i, 1] else u_id[i, 1] + u_id[i, 2] * t
    as.numeric(fixed_part + id_part)
  }

  cv_marker_i <- function(i, t) {
    # z_idm: (1, t) if Q_idm=2 else (1)
    # z_mk: (x1) if R_mk=1 else (empty)
    # Total CV_marker = mean over markers of [ mk_part + idm_part ]
    x1_val <- dataEvent$x1[i]
    lt <- length(t)
    per_marker_sum <- numeric(lt)

    for (d in seq_len(D)) {
      mk_part <- if (include_marker_only && R_mk > 0) x1_val * v_marker[d, 1] else 0
      idm_part <- if (Q_idm == 1) w_idscaled[i, d, 1] else w_idscaled[i, d, 1] + w_idscaled[i, d, 2] * t
      per_marker_sum <- per_marker_sum + mk_part + idm_part
    }
    per_marker_sum / D
  }
  cv_total_i <- function(i, t) cv_mean_i(i, t) + cv_marker_i(i, t)

  eta_w_i <- function(i) sum(haz_row(i) * gamma_w)

  hazard_i <- function(i, t) {
    h0(t) * exp(eta_w_i(i) + alpha_cv_total * cv_total_i(i, t))
  }

  cumhaz_i <- function(i, t) {
    if (t <= 0) {
      return(0)
    }
    out <- do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integration_control))
    as.numeric(out$value)
  }

  survival_prob_i <- function(i, t) exp(-cumhaz_i(i, t))

  draw_event_time <- function(i) {
    U <- runif(1)
    target_H <- -log(U)
    f <- function(t) cumhaz_i(i, t) - target_H
    br <- .find_bracket(f,
      lower = 0, upper = root_control$t_init, upper_max = root_control$t_max,
      expand = root_control$expand, max_expand = root_control$max_expand
    )
    if (is.null(br)) {
      return(list(time = time_cens, event = 0L, bracketing_failed = TRUE, U = U))
    }
    T <- uniroot(f, lower = br$lower, upper = br$upper)$root
    if (T > time_cens) {
      list(time = time_cens, event = 0L, bracketing_failed = FALSE, U = U)
    } else {
      list(time = T, event = 1L, bracketing_failed = FALSE, U = U)
    }
  }

  # Draw event times by inverse transform sampling
  death <- lapply(seq_len(n_id), draw_event_time)
  dataEvent$time <- vapply(death, `[[`, numeric(1), "time")
  dataEvent$event <- vapply(death, `[[`, integer(1), "event")

  # ---- longitudinal long format
  marker_levels <- paste0("m", seq_len(D))
  dataLong <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    do.call(rbind, lapply(seq_len(D), function(d) {
      t_max_i <- min(dataEvent$time[i], time_cens)
      t_obs <- sort(runif(n_t, 0, max(1e-8, t_max_i)))
      df <- data.frame(
        id = i,
        marker = factor(marker_levels[d], levels = marker_levels),
        time = t_obs,
        x1 = dataEvent$x1[i],
        x2 = dataEvent$x2[i],
        stringsAsFactors = FALSE
      )
      mu <- numeric(length(t_obs))
      for (k in seq_along(t_obs)) {
        tt <- t_obs[k]
        mu_fixed <- sum(fixed_row(i, tt) * beta)
        mu_id <- sum(id_row(tt) * u_id[i, ])
        mu_mk <- 0
        if (include_marker_only && R_mk > 0) mu_mk <- sum(mk_row(i, tt) * v_marker[d, ])
        mu_idm <- sum(idm_row(tt) * w_idscaled[i, d, ])
        mu[k] <- mu_fixed + mu_id + mu_mk + mu_idm
      }
      sig <- if (is.null(sigma_fun)) {
        rep(sigma_y, length(mu))
      } else {
        sigma_fun(i, t_obs, dataEvent)
      }
      if (is.null(sigma_fun) && !is.null(dist_formulas$sigma) && !is.null(beta_sigma)) {
        X_sigma <- .build_dist_matrix(dist_formulas$sigma, df, family_by_row = rep("student_t", nrow(df)))$X
        sig <- as.numeric(exp(X_sigma %*% beta_sigma))
      }
      if (length(sig) != length(mu)) {
        stop("sigma_fun must return a vector with length equal to times.")
      }
      nu_vec <- rep(nu, length(mu))
      if (!is.null(dist_formulas$nu) && !is.null(beta_nu)) {
        X_nu <- .build_dist_matrix(dist_formulas$nu, df, family_by_row = rep("student_t", nrow(df)))$X
        nu_vec <- 2 + exp(as.numeric(X_nu %*% beta_nu))
      }
      df$y <- vapply(seq_along(mu), function(k) mu[k] + rt(1, df = nu_vec[k]) * sig[k], numeric(1))
      df
    }))
  }))

  tmax <- max(dataEvent$time)

  # Store latent truth for debugging or benchmark comparisons
  truth <- list(
    beta = beta,
    gamma_w = gamma_w,
    alpha_cv_total = alpha_cv_total,
    nu = nu,
    sigma_y = sigma_y,
    beta_sigma = beta_sigma,
    beta_nu = beta_nu,
    formulaDist = dist_formulas,
    include_marker_only = include_marker_only,
    tau_u = tau_u,
    Lcorr_u = Lcorr_u,
    z_u = z_u,
    u_id = u_id,
    tau_v = tau_v,
    Lcorr_v = Lcorr_v,
    z_v = z_v,
    v_marker = v_marker,
    z_w_lat = z_w_lat,
    z_w = z_w,
    alpha_L = alpha_L,
    lambda_L = lambda_L,
    z_L = z_L,
    L_i = L_i,
    w_idscaled = w_idscaled
  )

  # Return helper closures for CV and hazard functions
  helpers <- list(
    cv_mean = cv_mean_i,
    cv_marker = cv_marker_i,
    cv_total = cv_total_i,
    hazard = hazard_i,
    cumhaz = cumhaz_i,
    survival_prob = survival_prob_i
  )

  list(dataLong = dataLong, dataEvent = dataEvent, truth = truth, helpers = helpers, tmax = tmax)
}

#' Simulate joint model data with formula-driven multistructure support
#'
#' @description
#' Generalised simulator for `joinme` that mirrors the fitting syntax as closely as
#' possible. The simulator supports:
#' - multivariate outcomes via marker-level families,
#' - multi-structure random effects from `formulaLong` (id, marker, marker-by-id),
#' - formula-based covariate generation through user-provided random generators,
#' - association-driven survival via `assoc`/`formulaAssoc` terms consistent with
#'   the model fit interface.
#'
#' The implementation is model-matrix based end-to-end, so every simulated component
#' is generated from the same formula machinery used during fitting.
#'
#' @param formulaLong Longitudinal formula (same role as in `joinme()`).
#'   Grouping terms may use `weighted(group, weights = <column>)` to mirror
#'   fitting syntax. The referenced weight column
#'   should be available in generated covariates (for example via
#'   `covariate_formulas`).
#' @param formulaEvent Event/survival formula (same role as in `joinme()`).
#' @param formulaVCov Optional covariance-regression formula for the id-specific
#'   marker-by-id covariance factor (same role as in `joinme_standata()`).
#'   This formula is evaluated on event-level covariates (one row per subject),
#'   must not include random-effect bars `( ... | ... )`, and must not include
#'   the longitudinal time variable.
#'
#'   Internally, this formula drives subject-specific lower-triangular entries
#'   of `L_i` used to scale marker-by-id latent effects; see `re_params$id_marker_cov`.
#' @param formulaDist Optional distributional regression formulas (same role as in `joinme()`).
#'   Supported LHS parameters are `sigma`, `nu`, `phi`, `alpha` (aliases:
#'   `alpha_skew`, `skew`), `phi_beta`, and `tau_sde`.
#'
#'   Family-scoped syntax is supported with square brackets:
#'   `param[family=<name>] ~ ...`, for example
#'   `sigma[family=student_t] ~ 1 + time`.
#'
#'   If no bracket is used (e.g. `sigma ~ 1 + time`), the formula applies to all
#'   rows where that parameter exists. If a bracket is used, it only applies to
#'   rows of that family and is estimated separately from other scopes.
#' @param formulaAssoc Optional formula for association terms in hazard, e.g. `~ cv_total + corr + vcov`.
#'   If provided, it overrides `assoc`.
#' @param transforms Optional transform specifications for association terms.
#'   Supports the same structure as `joinme(..., transforms=...)` for
#'   `cv_total`, `cv_mean`, `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.
#'   Simulation applies `cv_total` and `cv_marker` transforms at marker level before
#'   weighted averaging (aligned with fit/predict Stan semantics).
#'   Functional transforms support arithmetic and common nonlinear functions,
#'   including `inv_logit`/`expit`/`sigmoid`, `exp`, `log`, `sqrt`, `power`,
#'   `cbrt`, `softplus`/`log1p_exp`, trigonometric and hyperbolic functions.
#'
#'   Quick parameterisation reference:
#'   - `list(type = "identity")` or omitted term: identity transform (default).
#'   - `list(type = "functional", expr = ~ log1p(x))`: functional transform;
#'     `expr` is required.
#'   - `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
#'     direct monotone I-spline; `degree` defaults to `3` if omitted.
#'   - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
#'     penalised monotone I-spline; defaults are `n_knots = 6`, `degree = 3`,
#'     and `lambda = 1` when omitted.
#'   - `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
#'     monotone I-spline evaluated on `plogis(x)`; explicit knots are supplied
#'     on the expit scale.
#'   - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
#'     penalised monotone I-spline on `plogis(x)` in legacy plug-in mode.
#'   - `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`:
#'     piecewise-linear transform; `x` and `y` are required.
#'
#'   For monotone spline transforms during simulation:
#'   - `type = "ispline"`: provide `knots` and `coeff` directly (plus optional
#'     `degree`).
#'   - `type = "ispline_penalised"` (alias: `"ispline_penalized"`): provide
#'     training pairs `x` and `y`, plus
#'     smoothness penalty `lambda`. The simulator fits the monotone spline in R
#'     before evaluating the transformed association. The fitted plug-in spline
#'     uses anchored endpoint coefficients matching the Stan-estimated path:
#'     increasing splines run from `0` to `1`, while decreasing splines run from
#'     `1` to `0`. If `direction` is omitted, the simulator infers it from the
#'     supplied `(x, y)` pairs.
#'   - `type = "ispline_expit"` / `"ispline_expit_penalised"`: same semantics
#'     as the ordinary I-spline variants above, except the spline basis is built
#'     on `plogis(x)` for a bounded-domain representation. Training `x` values
#'     and explicit `knots` for these transform types are specified on that
#'     bounded expit scale. This is useful
#'     when the raw association feature has long tails or steep nonlinear effects.
#'
#'   In other words, simulation currently uses the legacy plug-in spline mode;
#'   it does not estimate spline coefficients jointly inside Stan.
#' @param marker_weights Optional base marker weights used as prior offsets for
#'   association aggregation. Effective weights are computed as
#'   `w_raw` and used for marker-averaged CV/CS terms.
#'   If `NULL`, equal weights are used and aggregated
#'   as weighted means divided by marker count. Inputs are treated as direct
#'   signed weights used in association aggregation.
#' @param n_id Number of subjects.
#' @param families Marker-specific family names.
#'   Use `jm_family()` entries to supply custom `link`/`inv_link` expressions.
#' @param marker_levels Optional marker names; defaults to `m1`, `m2`, ...
#' @param n_obs_per_marker_per_id Target number of observations per (id, marker).
#' @param times_obs Optional candidate observation time grid.
#' @param n_t Optional alias for `n_obs_per_marker_per_id` for compatibility.
#' @param seed RNG seed.
#' @param covariate_formulas Named or LHS formulas used to generate event-level covariates,
#'   e.g. `list(x1 ~ rnorm(n_id), x2 ~ rt(n_id, df = 5))`.
#' @param assoc Association components (same names as fit): `cv_total`, `cv_mean`,
#'   `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.
#'
#'   Meaning of each channel:
#'   - `cv_*`: current-value channels from longitudinal trajectories,
#'   - `cs_*`: current-slope channels from finite differences (`eps_cs`),
#'   - `corr`: off-diagonal correlation features from marker-by-id random effects.
#'   - `vcov`: lower-triangular Cholesky-factor entries from marker-by-id random
#'     effects, taken directly from the subject-specific `L` matrix.
#'
#'   You can also provide `formulaAssoc = ~ ...` to select channels; when present,
#'   it overrides `assoc`.
#' @param assoc_coefs Association coefficients for hazard terms.
#'   Non-`corr`/`vcov` terms accept scalar values. The `corr` term accepts a vector of
#'   off-diagonal correlation coefficients ordered as `(2,1), (3,1), (3,2), ...`
#'   in lower-triangular row-major order of the marker-by-id random-effect
#'   covariance dimension. The `vcov` term accepts lower-triangular Cholesky-factor
#'   entries ordered as `(1,1), (2,1), (2,2), ...)`; when `||` is used in the marker-by-id
#'   random-effects block, only diagonal `L` entries are used.
#'
#'   Accepted input forms:
#'   - named numeric vector, e.g. `c(cv_total = 0.4, cs_mean = -0.2)`,
#'   - named list, e.g. `list(cv_total = 0.4, corr = c(0.2, -0.1), vcov = c(0.4, 0.1, 0.5))`.
#'   Missing channels default to 0.
#' @param beta_long Fixed-effect coefficients for `formulaLong` fixed part. If NULL,
#'   coefficients are randomly generated and named by model-matrix columns.
#' @param beta_event Survival baseline-covariate coefficients for non-intercept
#'   terms in `formulaEvent` RHS.
#'   If NULL, coefficients are randomly generated.
#' @param dist_coefs Distributional fixed-effect coefficients for `formulaDist`
#'   parameters (`sigma`, `nu`, `phi`, `alpha`, `phi_beta`, `tau_sde`).
#'
#'   For each parameter, coefficients can be:
#'   - an unnamed numeric vector (matched by column order),
#'   - a named numeric vector (matched by model-matrix column names),
#'   - for family-scoped formulas, a named list with per-scope entries.
#'
#'   Family-scoped list syntax examples:
#'   - `dist_coefs = list(sigma = list(default = c("(Intercept)" = -0.3), gaussian = c(...), student_t = c(...)))`
#'   - alias keys like `"sigma[family='student_t']"` are also recognised and
#'     mapped to the matching family scope.
#' @param re_params Random-effects simulation controls.
#'
#'   Structure:
#'   - `id`: controls id-level random effects from `( ... | id)` in `formulaLong`.
#'   - `marker`: controls marker-level random effects from marker-only terms.
#'   - `id_marker_cov`: controls subject-specific covariance-regression for
#'     marker-by-id latent effects.
#'   - `dist`: controls random effects for distributional regressions in
#'     `formulaDist`.
#'
#'   For `id` and `marker`, each block is a list:
#'   - `sd`: scalar or length-K vector of random-effect standard deviations,
#'   - `corr`: KxK correlation matrix.
#'
#'   `id_marker_cov` fields:
#'   - `latent`: deprecated compatibility input. If supplied, its implied
#'     lower-triangular factor is folded into the baseline `alpha` intercepts
#'     before simulation. Marker-by-id latent seeds are still drawn as iid
#'     standard normal values. Prefer setting `alpha` directly in new code.
#'   - `alpha`: baseline linear predictors for entries of subject-specific lower
#'     triangular `L_i` (baseline when `formulaVCov` covariates and subject-level
#'     covariance perturbation are zero),
#'   - `beta`: covariate effects from `formulaVCov` design matrix (systematic
#'     subject-to-subject covariance shifts by observed covariates),
#'   - `lambda`: non-negative loading on an iid standard-normal subject latent
#'     perturbation; if a negative value is supplied, the simulator folds the
#'     sign into the latent draw so the effective model remains unchanged but
#'     follows the identified convention used during fitting,
#'   - `diag_link`: diagonal link for `L_i` diagonals (`"softplus"` or `"exp"`).
#'
#'   Element-wise covariance-regression form is:
#'   `eta_{im} = alpha_m + x_i^T beta_m + lambda_m z_{im}`, with
#'   `lambda_m >= 0` and `z_{im} ~ Normal(0, 1)`. Diagonal entries of `L_i` apply `diag_link`
#'   to keep them positive; off-diagonal entries remain on identity scale.
#'   Marker-by-id latent seeds are sampled as iid standard normal values and the
#'   final marker-by-id effects are obtained as `b_id = L_i z_id`.
#'   This means baseline and subject-specific covariance now both live in `L_i`,
#'   while `alpha`/`beta`/`lambda` control its level and heterogeneity.
#'
#'   Dimension rules for `id_marker_cov` entries follow marker-by-id random-effect
#'   dimension `Q_idm`:
#'   - if covariance is full: `M = Q_idm * (Q_idm + 1) / 2` lower-tri entries,
#'   - if covariance is forced diagonal (`||` in nested id-marker term): `M = Q_idm`.
#'
#'   `dist` block syntax:
#'   - `re_params$dist[[param]]` applies to all random-effect terms for that
#'     distributional parameter,
#'   - `re_params$dist[[param]]$terms[[j]]` optionally sets term-specific
#'     controls (same `sd`/`corr` fields as above).
#' @param family_params Family-specific simulation parameters.
#' @param h0 Optional baseline hazard function `h0(t)` for backward compatibility.
#'   If supplied, it takes precedence over `baseline_hazard`/`formulaBasehaz`.
#' @param baseline_hazard Optional baseline hazard specification. Supported forms:
#'   - character: one of `"constant"`, `"linear"`, `"piecewise"`, `"weibull"`, `"spline"`.
#'   - named list: `list(type = ..., ...)` with mode-specific parameters.
#' @param formulaBasehaz Optional formula-based baseline hazard model on time,
#'   e.g. `~ 1 + time + I(time^2)`.
#' @param beta_basehaz Optional coefficients for `formulaBasehaz` (aligned by
#'   model-matrix column names). If NULL, coefficients are generated.
#' @param time_cens Administrative censoring horizon.
#' @param eps_cs Finite-difference step for slope-type associations (`cs_*`).
#' @param integration_control Control list passed to `integrate()`.
#' @param root_control Root finding control for inverse-CDF sampling.
#' @param quadrature_nodes Optional Gauss-Kronrod node count for simulation-side
#'   survival integration. Allowed values: 7, 15, 31, 41, 51, 61.
#' @param n_workers Number of workers to use when `use_mirai = TRUE`.
#' @param use_mirai Logical; when TRUE and `n_workers > 1`, use `mirai` for
#'   parallel simulation if available. Parallel jobs are seeded deterministically
#'   from the main `seed` for reproducible simulations.
#' @param id_var,marker_var,time_var,y_var,event_time_var,event_var Column names aligned
#'   with `joinme_standata()` defaults.
#'
#' @return list(dataLong, dataEvent, true_params, marker_info, helpers, tmax)
#'
#' @examples
#' \dontrun{
#' sim <- simulate_joinme(
#'   n_id = 30,
#'   families = list(
#'     jm_family("gaussian"),
#'     jm_family("student_t", inv_link = ~ inv_logit(x / 2)),
#'     jm_family("skew_normal")
#'   ),
#'   quadrature_nodes = 31,
#'   formulaDist = list(
#'     sigma[family=gaussian] ~ 1 + x1,
#'     sigma[family=student_t] ~ 1 + time,
#'     nu[family=student_t] ~ 1,
#'     alpha_skew[family=skew_normal] ~ 1 + x1
#'   )
#' )
#' }
#' @export
simulate_joinme <- function(
  formulaLong = y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  formulaVCov = ~ 1,
  formulaDist = NULL,
  formulaAssoc = NULL,
  transforms = NULL,
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  n_obs_per_marker_per_id = 8,
  times_obs = seq(0, 5, length.out = 8),
  n_t = NULL,
  seed = 42,
  covariate_formulas = list(
    x1 ~ rnorm(n_id),
    x2 ~ rnorm(n_id)
  ),
  marker_weights = NULL,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6),
  beta_long = NULL,
  beta_event = NULL,
  dist_coefs = list(),
  re_params = list(
    id = list(sd = NULL, corr = NULL),
    marker = list(sd = NULL, corr = NULL),
    id_marker_cov = list(
      latent = list(sd = NULL, corr = NULL),
      alpha = NULL,
      beta = NULL,
      lambda = NULL,
      diag_link = "softplus"
    ),
    dist = list()
  ),
  family_params = list(
    gaussian = list(sigma = 1.0),
    student_t = list(sigma = 1.5, nu = 4),
    bernoulli = list(),
    binomial = list(trials = 10),
    poisson = list(),
    negbin2 = list(phi = 2),
    skew_normal = list(sigma = 1.0, alpha = 0),
    double_exponential = list(sigma = 1.0),
    skew_double_exponential = list(sigma = 1.0, tau_sde = 0.5),
    beta = list(phi_beta = 10),
    cumulative_logit = list(cutpoints = c(-1, 1))
  ),
  h0 = NULL,
  baseline_hazard = list(type = "weibull", shape = 1.4, scale = 6.0),
  formulaBasehaz = NULL,
  beta_basehaz = NULL,
  time_cens = 8.0,
  eps_cs = 1e-3,
  integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1.0, t_max = 50.0, expand = 1.7, max_expand = 60L),
  quadrature_nodes = 15L,
  n_workers = 1L,
  use_mirai = TRUE,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  y_var = "y",
  event_time_var = "time",
  event_var = "event"
) {
  # Workflow:
  # 1) Parse model formulas and random-effects structure.
  # 2) Generate subject-level covariates using formula-driven RNG expressions.
  # 3) Draw fixed/random coefficients and define latent trajectory evaluators.
  # 4) Simulate event times by inverse transform using hazard linked to assoc terms.
  # 5) Simulate observation times, generate responses by family, and return truth.
  set.seed(seed)
  assoc_coefs_missing <- missing(assoc_coefs)
  # Base seed for deterministic parallel jobs when mirai is enabled.
  seed_base <- as.integer(seed %||% 1L)
  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 1L
  use_mirai <- isTRUE(use_mirai) && n_workers > 1L
  if (use_mirai && !requireNamespace("mirai", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg mirai} not installed; falling back to serial simulation.")
    use_mirai <- FALSE
  }
  if (use_mirai) {
    mirai::daemons(n_workers)
    on.exit(mirai::daemons(0L), add = TRUE)
  }

  # Local bytecode evaluators used by simulation closures.
  # These are intentionally self-contained so mirai workers do not depend on
  # package namespace internals.
  .sim_eval_bytecode_scalar <- function(x, bytecode, const_data) {
    code <- as.integer(bytecode %||% integer(0))
    constants <- as.numeric(const_data %||% numeric(0))
    if (length(code) == 0L) return(as.numeric(x))

    stack <- numeric(0)
    const_idx <- 1L
    for (op in code) {
      if (op == 0L) {
        stack <- c(stack, as.numeric(x))
      } else if (op == 1L) {
        stack <- c(stack, constants[const_idx])
        const_idx <- const_idx + 1L
      } else if (op == 2L) {
        b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
        stack <- c(stack[-c(length(stack) - 1L, length(stack))], a + b)
      } else if (op == 3L) {
        b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
        stack <- c(stack[-c(length(stack) - 1L, length(stack))], a - b)
      } else if (op == 4L) {
        b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
        stack <- c(stack[-c(length(stack) - 1L, length(stack))], a * b)
      } else if (op == 5L) {
        b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
        stack <- c(stack[-c(length(stack) - 1L, length(stack))], a / b)
      } else if (op == 6L) {
        stack[length(stack)] <- log(stack[length(stack)])
      } else if (op == 7L) {
        stack[length(stack)] <- exp(stack[length(stack)])
      } else if (op == 8L) {
        stack[length(stack)] <- sqrt(stack[length(stack)])
      } else if (op == 9L) {
        stack[length(stack)] <- stats::plogis(stack[length(stack)])
      } else if (op == 10L) {
        stack[length(stack)] <- stats::qlogis(stack[length(stack)])
      } else if (op == 11L) {
        stack[length(stack)] <- 1 / stack[length(stack)]
      } else if (op == 12L) {
        b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
        stack <- c(stack[-c(length(stack) - 1L, length(stack))], a^b)
      } else if (op == 13L) {
        stack[length(stack)] <- sin(stack[length(stack)])
      } else if (op == 14L) {
        stack[length(stack)] <- cos(stack[length(stack)])
      } else if (op == 15L) {
        stack[length(stack)] <- tan(stack[length(stack)])
      } else if (op == 16L) {
        stack[length(stack)] <- abs(stack[length(stack)])
      } else if (op == 17L) {
        stack[length(stack)] <- stack[length(stack)]^2
      } else if (op == 18L) {
        stack[length(stack)] <- sinh(stack[length(stack)])
      } else if (op == 19L) {
        stack[length(stack)] <- cosh(stack[length(stack)])
      } else if (op == 20L) {
        stack[length(stack)] <- tanh(stack[length(stack)])
      } else if (op == 21L) {
        stack[length(stack)] <- asinh(stack[length(stack)])
      } else if (op == 22L) {
        stack[length(stack)] <- acosh(stack[length(stack)])
      } else if (op == 23L) {
        stack[length(stack)] <- atanh(stack[length(stack)])
      } else if (op == 24L) {
        stack[length(stack)] <- .softplus(stack[length(stack)])
      } else if (op == 25L) {
        a <- stack[length(stack)]
        stack[length(stack)] <- sign(a) * abs(a)^(1 / 3)
      } else if (op == 26L) {
        stack[length(stack)] <- stats::pnorm(stack[length(stack)])
      } else {
        cli::cli_abort("Unknown transform bytecode instruction: {op}.")
      }
    }
    stack[length(stack)]
  }

  .sim_eval_bytecode_vector <- function(x, bytecode, const_data) {
    code <- as.integer(bytecode %||% integer(0))
    constants <- as.numeric(const_data %||% numeric(0))
    vapply(as.numeric(x), .sim_eval_bytecode_scalar, numeric(1), bytecode = code, const_data = constants)
  }

  #' @keywords internal
  #' @param x Vector of inputs to process.
  #' @param fun Function to apply.
  #' @param ... Additional arguments passed to fun.
  #' @return List of results.
  .sim_parallel_lapply <- function(x, fun, ...) {
    # Parallel map with mirai; falls back to serial when mirai is disabled.
    # Each job sets a deterministic seed derived from the main seed.
    if (!use_mirai) return(lapply(x, fun, ...))
    args <- list(...)
    jobs <- lapply(seq_along(x), function(idx) {
      xi <- x[[idx]]
      seed_i <- seed_base + as.integer(idx)
      mirai::mirai({
        set.seed(seed_i)
        do.call(fun, c(list(xi), args))
      }, fun = fun, xi = xi, args = args, seed_i = seed_i)
    })
    # lapply(jobs, mirai::collect_mirai)
    mirai::collect_mirai(jobs, options = c('.stop', '.progress'))
  }

  .sim_align_coef <- function(col_names, user_coef = NULL, sd_default = 0.4, intercept_default = 0.0) {
    if (length(col_names) == 0) return(numeric(0))
    if (is.null(user_coef)) {
      out <- stats::rnorm(length(col_names), 0, sd_default)
      names(out) <- col_names
      if ("(Intercept)" %in% col_names) out["(Intercept)"] <- intercept_default
      return(out)
    }
    out <- rep(0, length(col_names))
    names(out) <- col_names
    if (!is.null(names(user_coef))) {
      keep <- intersect(names(user_coef), col_names)
      if (length(keep) > 0) {
        out[keep] <- as.numeric(user_coef[keep])
      }
    } else {
      out[seq_len(min(length(out), length(user_coef)))] <- as.numeric(user_coef)[seq_len(min(length(out), length(user_coef)))]
    }
    out
  }

  .sim_random_positive <- function(n, min_val = 0.25, max_val = 0.85) {
    stats::runif(n, min = min_val, max = max_val)
  }

  .sim_random_signed <- function(n, min_abs = 0.15, max_abs = 0.7) {
    signs <- sample(c(-1, 1), size = n, replace = TRUE)
    mags <- stats::runif(n, min = min_abs, max = max_abs)
    signs * mags
  }

  .sim_random_corr_matrix <- function(K, shrink_min = 0.15, shrink_max = 0.45) {
    if (K <= 0) return(matrix(0.0, 0, 0))
    if (K == 1) return(matrix(1.0, 1, 1))

    A <- matrix(stats::rnorm(K * K), nrow = K, ncol = K)
    C_rand <- stats::cov2cor(crossprod(A) + diag(K))
    shrink <- stats::runif(1, min = shrink_min, max = shrink_max)
    C <- (1 - shrink) * diag(K) + shrink * C_rand
    C <- 0.5 * (C + t(C))
    diag(C) <- 1.0
    C
  }

  .sim_resolve_re_block_cfg <- function(cfg, K, label) {
    if (K <= 0) {
      return(list(sd = numeric(0), corr = matrix(0.0, 0, 0), cov = matrix(0.0, 0, 0), Lcorr = matrix(0.0, 0, 0)))
    }

    sd_vec <- cfg$sd
    if (is.null(sd_vec)) {
      sd_vec <- .sim_random_positive(K)
    } else {
      sd_vec <- as.numeric(sd_vec)
      if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, K)
      if (length(sd_vec) != K || any(!is.finite(sd_vec)) || any(sd_vec <= 0)) {
        cli::cli_abort(c(
          x = "Random-effect SD specification is invalid for {label}.",
          i = "Expected {K} finite positive value(s)."
        ))
      }
    }

    corr <- cfg$corr
    if (is.null(corr)) {
      corr <- .sim_random_corr_matrix(K)
    } else {
      corr <- as.matrix(corr)
      if (!all(dim(corr) == c(K, K)) || any(!is.finite(corr))) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix is invalid for {label}.",
          i = "Expected a finite {K}x{K} matrix."
        ))
      }
      if (!isTRUE(all.equal(corr, t(corr), tolerance = 1e-8))) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must be symmetric for {label}.",
          i = "Ensure corr equals its transpose."
        ))
      }
      if (any(abs(diag(corr) - 1) > 1e-8)) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must have unit diagonal for {label}.",
          i = "All diagonal entries must equal 1."
        ))
      }
      if (any(abs(corr) > 1 + 1e-8)) {
        cli::cli_abort(c(
          x = "Random-effect correlation entries must lie in [-1, 1] for {label}.",
          i = "Check the supplied correlation matrix."
        ))
      }
      eig <- tryCatch(eigen(corr, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
      if (any(!is.finite(eig)) || min(eig) <= 1e-8) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must be positive definite for {label}.",
          i = "Supply a valid correlation matrix."
        ))
      }
    }

    cov_mat <- diag(as.numeric(sd_vec), K, K) %*% corr %*% diag(as.numeric(sd_vec), K, K)
    Lcorr <- t(chol(corr))

    list(sd = sd_vec, corr = corr, cov = cov_mat, Lcorr = Lcorr)
  }

  .sim_inverse_diag_link <- function(x, diag_link) {
    x <- as.numeric(x)
    if (identical(diag_link, "exp")) {
      return(log(x))
    }
    log(expm1(x))
  }

  .sim_strip_dist_prefix <- function(cols) {
    sub("^(all::|family=[^:]+::)", "", cols)
  }

  .sim_align_dist_coef <- function(col_names,
                                   user_coef = NULL,
                                   dist_scope = NULL,
                                   sd_default = 0.15,
                                   intercept_default = 0.0) {
    if (length(col_names) == 0) return(numeric(0))

    if (!is.list(user_coef) || is.null(dist_scope) || !.is_dist_scope(dist_scope)) {
      return(.sim_align_coef(col_names, user_coef, sd_default = sd_default, intercept_default = intercept_default))
    }

    out <- .sim_align_coef(col_names, NULL, sd_default = sd_default, intercept_default = intercept_default)
    default_keys <- c("default", "all", "allFamilies", "all_families")
    user_names <- names(user_coef) %||% character(0)

    # Optional default coefficients for all-scoped columns.
    key_default <- intersect(default_keys, user_names)
    if (length(key_default) > 0) {
      idx_all <- grepl("^all::", col_names)
      if (any(idx_all)) {
        cols_all <- .sim_strip_dist_prefix(col_names[idx_all])
        out[idx_all] <- .sim_align_coef(
          cols_all,
          user_coef[[key_default[1]]],
          sd_default = sd_default,
          intercept_default = intercept_default
        )
      }
    }

    # Family-specific coefficients; keys accepted: gaussian, student_t, family=gaussian.
    fam_keys <- setdiff(user_names, default_keys)
    for (k in fam_keys) {
      fam_raw <- sub("^family=", "", k)
      fam_name <- tryCatch(.canonical_family_name(fam_raw), error = function(e) fam_raw)
      idx_f <- grepl(paste0("^family=", fam_name, "::"), col_names)
      if (!any(idx_f)) next
      cols_f <- .sim_strip_dist_prefix(col_names[idx_f])
      out[idx_f] <- .sim_align_coef(
        cols_f,
        user_coef[[k]],
        sd_default = sd_default,
        intercept_default = intercept_default
      )
    }

    out
  }

  .sim_dist_coef_spec <- function(param_name, dist_coefs, dist_scope) {
    spec <- dist_coefs[[param_name]]
    if (is.null(dist_scope) || !.is_dist_scope(dist_scope)) return(spec)

    if (!is.list(spec) || is.null(names(spec))) {
      out <- if (is.null(spec)) list() else list(default = spec)
    } else {
      out <- spec
    }

    fam_names <- names(dist_scope$by_family %||% list())
    if (length(fam_names) == 0) return(out)

    for (fam in fam_names) {
      pat <- paste0("^", param_name, "\\[family=['\"]?", fam, "['\"]?\\]$")
      hit <- grep(pat, names(dist_coefs), perl = TRUE)
      if (length(hit) > 0) {
        out[[fam]] <- dist_coefs[[names(dist_coefs)[hit[1]]]]
      }
    }
    out
  }

  .sim_eval_covariate_formulas <- function(n_subjects, formulas, base_df) {
    if (is.null(formulas) || length(formulas) == 0) return(base_df)
    if (!is.list(formulas)) formulas <- list(formulas)
    if (is.null(names(formulas))) names(formulas) <- rep("", length(formulas))

    out <- base_df
    for (j in seq_along(formulas)) {
      fj <- formulas[[j]]
      if (!inherits(fj, "formula")) fj <- stats::as.formula(fj)

      lhs_name <- names(formulas)[j]
      if (!nzchar(lhs_name)) {
        if (length(fj) < 3) {
          cli::cli_abort(c(
            x = "Each entry in {.arg covariate_formulas} must be named or have an LHS.",
            i = "Example: list(x1 ~ rnorm(n_id), x2 ~ runif(n_id))."
          ))
        }
        lhs_vars <- all.vars(fj[[2]])
        if (length(lhs_vars) != 1) {
          cli::cli_abort(c(
            x = "Covariate generation LHS must contain exactly one variable.",
            i = "Use syntax like x1 ~ rnorm(n_id)."
          ))
        }
        lhs_name <- lhs_vars[1]
      }

      rhs_expr <- if (length(fj) >= 3) fj[[3]] else fj[[2]]
      eval_env <- list2env(c(as.list(out), list(n_id = n_subjects, id = seq_len(n_subjects))), parent = parent.frame())
      sampled <- eval(rhs_expr, envir = eval_env)

      if (length(sampled) == 1) sampled <- rep(sampled, n_subjects)
      if (length(sampled) != n_subjects) {
        cli::cli_abort(c(
          x = "Covariate formula for '{lhs_name}' returned length {length(sampled)}.",
          i = "Expected length {.val {n_subjects}} or scalar."
        ))
      }
      out[[lhs_name]] <- sampled
    }
    out
  }

  .sim_draw_re_block <- function(n_group, K, cfg, label, resolved = NULL) {
    if (K <= 0) return(matrix(0.0, n_group, 0))
    resolved <- resolved %||% .sim_resolve_re_block_cfg(cfg, K, label)
    sd_vec <- resolved$sd
    corr <- resolved$corr
    Sigma <- diag(as.numeric(sd_vec), K, K) %*% corr %*% diag(as.numeric(sd_vec), K, K)
    R <- chol(Sigma)
    Z <- matrix(stats::rnorm(n_group * K), n_group, K)
    Z %*% R
  }

  .sim_rhs_matrix <- function(rhs_list, data) {
    if (length(rhs_list) == 0) {
      return(matrix(0.0, nrow(data), 0))
    }
    mats <- lapply(rhs_list, function(rhs) .mm(rhs, data))
    do.call(cbind, mats)
  }

  .sim_detect_time_cols <- function(mat_builder, prototype_df, time_var, tol = 1e-10) {
    mat_ref <- mat_builder(prototype_df)
    if (is.null(mat_ref) || ncol(mat_ref) == 0L) return(integer(0))

    probe_df <- prototype_df
    probe_df[[time_var]] <- probe_df[[time_var]] + 1
    mat_probe <- mat_builder(probe_df)
    changed <- colSums(abs(mat_probe - mat_ref) > tol) > 0
    which(changed)
  }

  .sim_get_family_param <- function(fam_name, param_name, fallback) {
    val <- family_params[[fam_name]][[param_name]]
    if (is.null(val)) fallback else val
  }

  .sim_random_assoc_parts <- function(assoc, M_corr, M_vcov) {
    scalar <- assoc_coef_scalar
    scalar_terms <- intersect(names(scalar), assoc)
    if (length(scalar_terms) > 0) {
      for (term in scalar_terms) {
        scalar[[term]] <- if (term %in% c("cv_total", "cv_marker", "cs_total", "cs_marker")) {
          stats::runif(1, min = 0.15, max = 0.75)
        } else {
          .sim_random_signed(1, min_abs = 0.1, max_abs = 0.75)
        }
      }
    }
    list(
      scalar = scalar,
      corr = if (M_corr > 0 && "corr" %in% assoc) .sim_random_signed(M_corr, min_abs = 0.1, max_abs = 0.5) else rep(0.0, M_corr),
      vcov = if (M_vcov > 0 && "vcov" %in% assoc) .sim_random_signed(M_vcov, min_abs = 0.1, max_abs = 0.8) else rep(0.0, M_vcov)
    )
  }

  #' @keywords internal
  #' @param x Numeric vector on link scale.
  #' @param inv_link_bc Parsed inverse-link bytecode.
  #' @return Numeric vector on response scale.
  .sim_apply_inv_link_bc <- function(x, inv_link_bc) {
    # Apply a parsed inverse-link bytecode using the shared interpreter.
    .sim_eval_bytecode_vector(
      x = x,
      bytecode = inv_link_bc$bytecode %||% inv_link_bc$opcodes,
      const_data = inv_link_bc$const_data %||% numeric(0)
    )
  }

  #' @keywords internal
  #' @param families Family specification input (vector or list).
  #' @param D Number of markers.
  #' @return List with family codes, link names, and inverse-link bytecode per marker.
  .sim_parse_family_specs <- function(families, D) {
    # Normalise family/link specs to a per-marker list of bytecode maps.
    if (length(families) == 1L && !is.list(families)) {
      families <- rep(list(families), D)
    } else if (!is.list(families)) {
      families <- as.list(families)
    }
    if (length(families) != D) {
      cli::cli_abort("families must have length {D} (number of markers), got {length(families)}")
    }

    specs <- lapply(families, .extract_family_and_link)
    family_codes <- vapply(specs, function(s) s$family_code, integer(1))
    link_names <- vapply(specs, function(s) s$link_name, character(1))
    inv_link_specs <- lapply(specs, function(s) s$inv_link_bc)

    list(
      family_codes = as.integer(family_codes),
      link_names = as.character(link_names),
      inv_link_specs = inv_link_specs
    )
  }

  .sim_sample_ordinal <- function(eta, cutpoints) {
    cp <- sort(as.numeric(cutpoints))
    cdf_vals <- stats::plogis(cp - eta)
    probs <- c(cdf_vals[1], diff(cdf_vals), 1 - cdf_vals[length(cdf_vals)])
    probs <- pmax(probs, 1e-12)
    probs <- probs / sum(probs)
    sample.int(length(probs), size = 1, prob = probs)
  }

  .sim_assoc_from_formula <- function(f) {
    if (is.null(f)) return(NULL)
    vars <- all.vars(f)
    unique(vars[vars %in% c("cv_total", "cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker", "corr", "vcov")])
  }

  .sim_normalize_transform_list <- function(transform_list) {
    unclass(.normalise_joinme_tf_input(transform_list, validate = FALSE))
  }

  .sim_canonicalise_transform_type <- function(type) {
    type <- as.character(type %||% "identity")[1]
    if (identical(type, "ispline_penalized")) {
      return("ispline_penalised")
    }
    if (identical(type, "ispline_expit_penalized")) {
      return("ispline_expit_penalised")
    }
    if (identical(type, "ispline_exp_penalised")) {
      return("ispline_expit_penalised")
    }
    if (identical(type, "ispline_exp_penalized")) {
      return("ispline_expit_penalised")
    }
    type
  }

  .sim_is_ispline_transform_type <- function(spec_or_type) {
    type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
    .sim_canonicalise_transform_type(type) %in% c(
      "ispline",
      "ispline_penalised",
      "pmonospline",
      "pmono",
      "ispline_expit",
      "ispline_expit_penalised"
    )
  }

  .sim_transform_uses_expit_input <- function(spec_or_type) {
    type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
    .sim_canonicalise_transform_type(type) %in% c("ispline_expit", "ispline_expit_penalised")
  }

  .sim_validate_expit_domain_values <- function(values, arg_name) {
    values <- as.numeric(values)
    if (!length(values)) {
      return(values)
    }
    if (any(!is.finite(values))) {
      cli::cli_abort(c(
        x = "{.arg {arg_name}} must contain only finite numeric values.",
        i = "For {.val ispline_expit} and {.val ispline_expit_penalised}, provide spline inputs on the expit scale in [0, 1]."
      ))
    }
    if (any(values < 0 | values > 1)) {
      cli::cli_abort(c(
        x = "{.arg {arg_name}} must lie on the expit scale for expit-based spline transforms.",
        i = "Supply values in [0, 1]; {.fn simulate_joinme} applies {.fn plogis} only to the raw association feature being transformed."
      ))
    }
    values
  }

  .sim_transform_input_for_spec <- function(x, spec_or_type) {
    x <- as.numeric(x)
    if (.sim_transform_uses_expit_input(spec_or_type)) {
      return(stats::plogis(x))
    }
    x
  }

  .sim_resolve_monotone_direction <- function(direction = NULL, default = "increasing") {
    direction <- direction %||% default
    if (is.numeric(direction) && length(direction) == 1L && is.finite(direction)) {
      if (direction > 0) return(1L)
      if (direction < 0) return(-1L)
    }
    direction_chr <- tolower(trimws(as.character(direction)[1]))
    if (direction_chr %in% c("increasing", "increase", "inc", "+1", "1")) return(1L)
    if (direction_chr %in% c("decreasing", "decrease", "dec", "-1")) return(-1L)
    cli::cli_abort(c(
      x = "{.arg direction} must be either {.val increasing} or {.val decreasing}.",
      i = "Numeric aliases {.val 1} and {.val -1} are also accepted."
    ))
  }

  .sim_infer_monotone_direction <- function(x, y) {
    ord <- order(x)
    x <- as.numeric(x)[ord]
    y <- as.numeric(y)[ord]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(y) < 2L) return(1L)
    if (y[length(y)] < y[1]) return(-1L)
    1L
  }

  .sim_make_penalised_ispline_transform <- function(spec) {
    if (is.null(spec$x) || is.null(spec$y)) {
      cli::cli_abort(c(
        x = "Penalised I-spline requires {.arg x} and {.arg y}.",
        i = "Provide input-output pairs for the monotone transform."
      ))
    }

    x <- as.numeric(spec$x)
    if (.sim_transform_uses_expit_input(spec)) {
      x <- .sim_validate_expit_domain_values(x, "x")
    }
    y <- as.numeric(spec$y)
    if (length(x) != length(y) || length(x) < 2L) {
      cli::cli_abort(c(
        x = "{.arg x} and {.arg y} must have the same length >= 2.",
        i = "Check the transform training data."
      ))
    }

    spline_direction <- .sim_resolve_monotone_direction(spec$direction, default = .sim_infer_monotone_direction(x, y))

    lambda <- as.numeric(spec$lambda %||% 1.0)
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
      cli::cli_abort(c(
        x = "{.arg lambda} must be a non-negative numeric scalar.",
        i = "Example: lambda = 1.0."
      ))
    }

    degree <- as.integer(spec$degree %||% 3L)
    if (degree < 1L) {
      cli::cli_abort(c(
        x = "{.arg degree} must be >= 1.",
        i = "Typical choice: degree = 3 (cubic)."
      ))
    }

    knots_raw <- spec$knots
    if (is.null(knots_raw)) {
      n_knots <- as.integer(spec$n_knots %||% 6L)
      if (n_knots < 2L) {
        cli::cli_abort(c(
          x = "{.arg n_knots} must be >= 2.",
          i = "Include boundary knots at the ends."
        ))
      }
      probs <- seq(0, 1, length.out = n_knots)
      knots <- as.numeric(stats::quantile(x, probs = probs, names = FALSE))
    } else {
      knots <- as.numeric(knots_raw)
      if (.sim_transform_uses_expit_input(spec)) {
        knots <- .sim_validate_expit_domain_values(knots, "knots")
      }
    }

    if (length(knots) < 2L || any(diff(knots) <= 0)) {
      cli::cli_abort(c(
        x = "{.arg knots} must be a strictly increasing vector with at least 2 values.",
        i = "Include boundary knots at the ends."
      ))
    }

    if (!requireNamespace("splines2", quietly = TRUE)) {
      cli::cli_abort(c(
        x = "Package {.pkg splines2} is required for penalised I-splines.",
        i = "Install splines2 or provide explicit coeff/knots for type = 'ispline'."
      ))
    }

    internal_knots <- if (length(knots) > 2L) knots[2:(length(knots) - 1L)] else numeric(0)
    boundary_knots <- c(knots[1], knots[length(knots)])
    basis <- splines2::iSpline(
      x,
      knots = internal_knots,
      degree = degree,
      intercept = TRUE,
      Boundary.knots = boundary_knots
    )
    basis <- as.matrix(basis)

    weights <- spec$weights
    if (!is.null(weights)) {
      weights <- as.numeric(weights)
      if (length(weights) != length(x) || any(weights < 0)) {
        cli::cli_abort(c(
          x = "{.arg weights} must be non-negative and the same length as {.arg x}.",
          i = "Remove weights or provide a valid vector."
        ))
      }
    } else {
      weights <- rep(1, length(x))
    }

    diff_order <- as.integer(spec$diff_order %||% 2L)
    if (diff_order < 1L) {
      cli::cli_abort(c(
        x = "{.arg diff_order} must be >= 1.",
        i = "Typical choice: diff_order = 2."
      ))
    }

    n_coef <- ncol(basis)
    if (n_coef < 2L) {
      cli::cli_abort(c(
        x = "Penalised I-spline basis must have at least 2 coefficients.",
        i = "Increase the number of knots or the spline degree."
      ))
    }
    Dmat <- diff(diag(n_coef), differences = diff_order)

    w_sqrt <- sqrt(weights)
    B_w <- basis * w_sqrt
    y_w <- y * w_sqrt

    init_raw <- tryCatch(
      as.numeric(qr.solve(crossprod(B_w), crossprod(B_w, y_w))),
      error = function(e) rep(0, n_coef)
    )

    .sim_anchored_coeff_from_z <- function(z) {
      z <- as.numeric(z)
      z <- z - max(z)
      delta <- exp(z)
      delta <- delta / sum(delta)
      c(0, cumsum(delta))
    }

    init_coeff <- pmax(0, init_raw)
    init_coeff <- init_coeff - init_coeff[1]
    span <- max(init_coeff, na.rm = TRUE)
    if (!is.finite(span) || span <= 0) {
      init_coeff <- seq(0, 1, length.out = n_coef)
    } else {
      init_coeff <- init_coeff / span
      init_coeff <- pmin(pmax(init_coeff, 0), 1)
      init_coeff <- cummax(init_coeff)
      init_coeff[1] <- 0
      init_coeff[n_coef] <- 1
    }
    delta_init <- diff(init_coeff)
    delta_init <- pmax(delta_init, .Machine$double.eps)
    delta_init <- delta_init / sum(delta_init)
    z_init <- log(delta_init)

    fn <- function(z) {
      b <- .sim_anchored_coeff_from_z(z)
      fit_vals <- as.numeric(B_w %*% b)
      if (spline_direction < 0L) {
        fit_vals <- w_sqrt * (sum(b) - as.numeric(basis %*% b))
      }
      r <- y_w - fit_vals
      pen <- if (lambda > 0 && nrow(Dmat) > 0) Dmat %*% b else 0
      0.5 * sum(r^2) + 0.5 * lambda * sum(pen^2)
    }

    opt <- optim(z_init, fn, method = "BFGS", control = list(maxit = 1000))
    if (opt$convergence != 0) {
      cli::cli_warn(c(
        x = "Penalised I-spline optimisation did not fully converge.",
        i = "Consider increasing lambda or adjusting knots."
      ))
    }

    list(
      type = if (.sim_transform_uses_expit_input(spec)) "ispline_expit" else "ispline",
      knots = knots,
      raw_knots = if (is.null(knots_raw)) NULL else as.numeric(knots_raw),
      coeff = as.numeric(.sim_anchored_coeff_from_z(opt$par)),
      degree = degree,
      direction = if (spline_direction < 0L) "decreasing" else "increasing",
      spline_direction = spline_direction
    )
  }

  .sim_make_assoc_transform <- function(spec, term_name) {
    if (is.null(spec) || is.null(spec$type) || identical(spec$type, "identity")) {
      return(function(x) x)
    }

    tf_type <- .sim_canonicalise_transform_type(spec$type)

    if (identical(tf_type, "functional")) {
      bc <- parse_transform_expr(spec$expr)
      return(function(x) {
        .sim_eval_bytecode_vector(
          x = x,
          bytecode = bc$bytecode %||% bc$opcodes,
          const_data = bc$const_data %||% numeric(0)
        )
      })
    }

    if (.sim_is_ispline_transform_type(tf_type)) {
      if (tf_type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised")) {
        if (is.null(spec$y)) {
          cli::cli_abort(c(
            x = "Simulation currently requires both {.arg x} and {.arg y} for penalised spline transforms.",
            i = "Use legacy plug-in mode in {.fn simulate_joinme} by supplying x/y pairs, or provide an explicit {.val ispline} transform instead.",
            i = "Stan-estimated penalised splines without y are supported in {.fn joinme}, not in simulation."
          ))
        }
        spec <- .sim_make_penalised_ispline_transform(spec)
      }
      if (!requireNamespace("splines2", quietly = TRUE)) {
        cli::cli_abort(c(
          x = "Transform type {.val {tf_type}} for {.val {term_name}} requires {.pkg splines2}.",
          i = "Install {.pkg splines2} or use identity/functional/pwlin in simulation."
        ))
      }
      knots <- as.numeric(spec$knots)
      if (.sim_transform_uses_expit_input(spec)) {
        knots <- .sim_validate_expit_domain_values(knots, "knots")
      }
      coeff <- as.numeric(spec$coeff)
      degree <- as.integer(spec$degree %||% 3L)
      return(function(x) {
        x <- .sim_transform_input_for_spec(x, spec)
        if (length(knots) < 2) {
          cli::cli_abort(c(
            x = "I-spline transform for {.val {term_name}} requires at least 2 knots.",
            i = "Provide boundary knots as the first and last entries."
          ))
        }
        boundary <- c(knots[1], knots[length(knots)])
        x <- pmin(pmax(x, boundary[1]), boundary[2])
        internal_knots <- if (length(knots) > 2) knots[2:(length(knots) - 1)] else numeric(0)
        basis <- splines2::iSpline(x,
          knots = internal_knots,
          degree = degree,
          intercept = TRUE,
          Boundary.knots = boundary
        )
        coef_len <- ncol(basis)
        if (length(coeff) < coef_len) {
          coeff_use <- c(coeff, rep(0, coef_len - length(coeff)))
        } else {
          coeff_use <- coeff[seq_len(coef_len)]
        }
        vals <- as.numeric(basis %*% coeff_use)
        if (.sim_resolve_monotone_direction(spec$direction %||% spec$spline_direction %||% 1L, default = 1L) < 0L) {
          vals <- sum(coeff_use) - vals
        }
        vals
      })
    }

    if (identical(tf_type, "pwlin")) {
      xk <- as.numeric(spec$x)
      yk <- as.numeric(spec$y)
      if (length(xk) < 2 || length(yk) != length(xk)) {
        cli::cli_abort(c(
          x = "Piecewise transform for {.val {term_name}} needs matching x/y vectors (length >= 2).",
          i = "Use list(type='pwlin', x=..., y=...)."
        ))
      }
      return(function(x) {
        stats::approx(x = xk, y = yk, xout = as.numeric(x), method = "linear", rule = 2)$y
      })
    }

    cli::cli_abort(c(
      x = "Unsupported transform type for {.val {term_name}}: {.val {tf_type}}.",
      i = "Use one of identity, functional, ispline, ispline_penalised, ispline_expit, ispline_expit_penalised, pwlin."
    ))
  }

  .sim_make_baseline_hazard <- function(spec, formula_basehaz = NULL, beta_basehaz = NULL) {
    if (!is.null(formula_basehaz)) {
      f_bh <- stats::as.formula(formula_basehaz)
      return(function(t) {
        df_t <- data.frame(time = as.numeric(t))
        X_bh <- stats::model.matrix(f_bh, data = df_t)
        coef_bh <- .sim_align_coef(
          colnames(X_bh),
          user_coef = beta_basehaz,
          sd_default = 0.2,
          intercept_default = -2.0
        )
        as.numeric(exp(X_bh %*% coef_bh))
      })
    }

    if (is.null(spec)) {
      spec <- list(type = "weibull", shape = 1.4, scale = 6.0)
    }
    if (is.character(spec)) {
      spec <- list(type = spec)
    }
    if (!is.list(spec) || is.null(spec$type)) {
      cli::cli_abort(c(
        x = "{.arg baseline_hazard} must be a character mode or named list with {.arg type}.",
        i = "Supported types: constant, linear, piecewise, weibull, spline."
      ))
    }

    mode <- tolower(as.character(spec$type)[1])
    if (mode == "constant") {
      rate <- as.numeric(spec$rate %||% spec$lambda %||% 0.08)
      return(function(t) rep(rate, length(t)))
    }

    if (mode == "linear") {
      intercept <- as.numeric(spec$intercept %||% -2.2)
      slope <- as.numeric(spec$slope %||% 0.25)
      return(function(t) {
        t <- as.numeric(t)
        .softplus(intercept + slope * pmax(t, 0))
      })
    }

    if (mode %in% c("piecewise", "pwlin", "piecewise_linear")) {
      breaks <- sort(as.numeric(spec$breaks %||% c(2, 5)))
      rates <- as.numeric(spec$rates %||% c(0.05, 0.12, 0.25))
      if (length(rates) != length(breaks) + 1L) {
        cli::cli_abort(c(
          x = "Piecewise baseline requires length(rates) = length(breaks) + 1.",
          i = "Got {.val {length(rates)}} rates and {.val {length(breaks)}} breaks."
        ))
      }
      return(function(t) {
        t <- pmax(as.numeric(t), 0)
        idx <- findInterval(t, vec = breaks, rightmost.closed = TRUE) + 1L
        rates[idx]
      })
    }

    if (mode == "weibull") {
      shape <- as.numeric(spec$shape %||% 1.4)
      scale <- as.numeric(spec$scale %||% 6.0)
      return(weibull_h0(shape = shape, scale = scale))
    }

    if (mode %in% c("spline", "bs", "ns")) {
      basis_type <- tolower(as.character(spec$basis %||% if (mode == "ns") "ns" else "bs")[1])
      knots <- as.numeric(spec$knots %||% stats::quantile(seq(0, time_cens, length.out = 100), probs = c(0.25, 0.5, 0.75)))
      degree <- as.integer(spec$degree %||% 3L)
      coef <- spec$coef
      intercept <- as.logical(spec$intercept %||% TRUE)

      return(function(t) {
        t <- pmax(as.numeric(t), 0)
        if (basis_type == "ns") {
          B <- splines::ns(t, knots = knots, Boundary.knots = c(0, time_cens), intercept = intercept)
        } else {
          B <- splines::bs(t, knots = knots, Boundary.knots = c(0, time_cens), degree = degree, intercept = intercept)
        }
        B <- as.matrix(B)
        coef_vec <- .sim_align_coef(
          colnames(B),
          user_coef = coef,
          sd_default = 0.25,
          intercept_default = -2.2
        )
        as.numeric(exp(B %*% coef_vec))
      })
    }

    cli::cli_abort(c(
      x = "Unsupported baseline hazard mode: {.val {mode}}.",
      i = "Use one of: constant, linear, piecewise, weibull, spline."
    ))
  }

  # ---- Basic dimensions and marker metadata
  if (!is.null(n_t)) n_obs_per_marker_per_id <- n_t
  if (length(families) == 1L && !is.null(marker_levels) && length(marker_levels) > 1L) {
    families <- rep(families, length(marker_levels))
  }
  D <- length(families)
  if (D < 1) {
    cli::cli_abort(c(
      x = "{.arg families} must include at least one marker family.",
      i = "Example: families = c('gaussian', 'student_t')."
    ))
  }
  if (is.null(marker_levels)) marker_levels <- paste0("m", seq_len(D))
  if (length(marker_levels) != D) {
    cli::cli_abort(c(
      x = "{.arg marker_levels} length must match {.arg families}.",
      i = "Expected {D}, got {length(marker_levels)}."
    ))
  }
  # ---- Resolve family specs and inverse-link bytecode per marker
  family_spec <- .sim_parse_family_specs(families, D)
  family_codes <- family_spec$family_codes
  link_names <- family_spec$link_names
  inv_link_specs <- family_spec$inv_link_specs
  family_names <- vapply(family_codes, .family_code_to_name, character(1))

  marker_weights_raw <- marker_weights %||% rep(1, D)
  marker_weights_raw <- as.numeric(unlist(marker_weights_raw, use.names = FALSE))
  if (length(marker_weights_raw) == 1L) marker_weights_raw <- rep(marker_weights_raw, D)
  if (length(marker_weights_raw) != D) {
    if (length(marker_weights_raw) > D) {
      marker_weights_raw <- marker_weights_raw[seq_len(D)]
    } else {
      marker_weights_raw <- rep(marker_weights_raw, length.out = D)
    }
  }
  if (any(!is.finite(marker_weights_raw))) {
    cli::cli_abort(c(
      x = "{.arg marker_weights} must be finite numeric with length equal to number of markers ({D}).",
      i = "Provide one weight per marker or a scalar recycled across markers."
    ))
  }
  marker_weights_eff <- marker_weights_raw

  # ---- Parse longitudinal formula structure
  f_exp <- reformulas::expandDoubleVerts(formulaLong)
  bars <- reformulas::findbars(f_exp)
  if (length(bars) == 0) {
    cli::cli_abort(c(
      x = "{.arg formulaLong} must include random-effects terms.",
      i = "Include at least ( ... | id ) and typically a marker block."
    ))
  }
  fixed_formula <- reformulas::nobars(f_exp)
  fixed_rhs <- stats::update(fixed_formula, . ~ .)
  fixed_rhs[[2]] <- NULL

  grp_names <- vapply(bars, function(b) .group_name_from_expr(b[[3]]), character(1))
  indep_flags <- .resolve_re_independence(formulaLong, marker_var = marker_var, id_var = id_var)
  id_idx <- which(grp_names == id_var)
  id_rhs_list <- if (length(id_idx) > 0) .bar_terms_to_rhs_list(bars[id_idx]) else list()

  nested_terms <- .extract_nested_marker_terms(formulaLong, marker_var = marker_var, id_var = id_var)
  mk_rhs_list <- nested_terms$mk_rhs_list
  idm_rhs_list <- nested_terms$idm_rhs_list

  # ---- Build event-level frame and covariates
  dataEvent <- data.frame(id = seq_len(n_id), stringsAsFactors = FALSE)
  dataEvent <- .sim_eval_covariate_formulas(n_subjects = n_id, formulas = covariate_formulas, base_df = dataEvent)
  names(dataEvent)[names(dataEvent) == "id"] <- id_var

  # Covariance-regression design for the subject-specific marker-by-id covariance.
  formulaVCov <- .resolve_vcov_formula(
    formulaVCov = formulaVCov,
    default = ~ 1,
    context = "simulate_joinme()"
  )
  vcov_design <- .build_vcov_design(
    formulaVCov = formulaVCov,
    dataEvent = dataEvent,
    time_var = time_var,
    context = "simulate_joinme()"
  )
  K_cov <- vcov_design$K_cov
  Xcov <- vcov_design$Xcov

  # ---- Build prototype data for matrix column alignment
  prototype <- dataEvent[rep(1, D), , drop = FALSE]
  prototype[[marker_var]] <- factor(marker_levels, levels = marker_levels)
  prototype[[time_var]] <- rep(mean(range(times_obs)), D)

  X_proto <- .mm(fixed_rhs, prototype)
  Z_id_proto <- .sim_rhs_matrix(id_rhs_list, prototype)
  Z_mk_proto <- .sim_rhs_matrix(mk_rhs_list, prototype)
  Z_idm_proto <- .sim_rhs_matrix(idm_rhs_list, prototype)

  idx_time_beta <- .sim_detect_time_cols(function(df) .mm(fixed_rhs, df), prototype, time_var)
  idx_time_uid <- .sim_detect_time_cols(function(df) .sim_rhs_matrix(id_rhs_list, df), prototype, time_var)
  idx_time_vmk <- .sim_detect_time_cols(function(df) .sim_rhs_matrix(mk_rhs_list, df), prototype, time_var)
  idx_time_widm <- .sim_detect_time_cols(function(df) .sim_rhs_matrix(idm_rhs_list, df), prototype, time_var)

  beta_long <- .sim_align_coef(colnames(X_proto), beta_long, sd_default = 0.35, intercept_default = 1.0)

  assoc_from_formula <- .sim_assoc_from_formula(formulaAssoc)
  assoc_effective <- if (!is.null(assoc_from_formula)) assoc_from_formula else assoc
  assoc_effective <- unique(assoc_effective)
  if (length(assoc_effective) == 0) assoc_effective <- "cv_total"

  W_event <- .mm_event(formulaEvent, dataEvent)

  if (is.null(beta_event)) {
    beta_event <- .sim_align_coef(
      colnames(W_event),
      user_coef = NULL,
      sd_default = 0.25,
      intercept_default = 0.0
    )
  } else {
    beta_event <- .sim_align_coef(
      colnames(W_event),
      beta_event,
      sd_default = 0.25,
      intercept_default = 0.0
    )
  }

  # ---- Draw random effects from user-configurable covariance structures
  K_id <- ncol(Z_id_proto)
  K_mk <- ncol(Z_mk_proto)
  K_idm <- ncol(Z_idm_proto)

  re_id_effective <- .sim_resolve_re_block_cfg(re_params$id %||% list(), K_id, "id")
  re_marker_effective <- .sim_resolve_re_block_cfg(re_params$marker %||% list(), K_mk, "marker")
  re_id <- .sim_draw_re_block(n_id, K_id, re_params$id %||% list(), "id", resolved = re_id_effective)
  re_marker <- .sim_draw_re_block(D, K_mk, re_params$marker %||% list(), "marker", resolved = re_marker_effective)
  re_cov_cfg <- re_params[["id_marker_cov", exact = TRUE]] %||% list()
  re_idm_cfg <- re_cov_cfg[["latent", exact = TRUE]]
  re_idm_legacy <- NULL
  latent_compat_factor <- if (K_idm > 0) diag(K_idm) else matrix(0.0, 0, 0)
  if (!is.null(re_idm_cfg) && K_idm > 0) {
    re_idm_legacy <- .sim_resolve_re_block_cfg(re_idm_cfg, K_idm, "id_marker_cov$latent")
    latent_sigma <- diag(re_idm_legacy$sd, K_idm, K_idm) %*% re_idm_legacy$corr %*% diag(re_idm_legacy$sd, K_idm, K_idm)
    latent_compat_factor <- t(chol(latent_sigma))
  }
  re_idm_effective <- list(
    mode = "iid_standard_normal",
    sd = if (K_idm > 0) rep(1, K_idm) else numeric(0),
    corr = if (K_idm > 0) diag(K_idm) else matrix(0.0, 0, 0),
    legacy_input = re_idm_legacy
  )
  re_idm_flat <- if (K_idm > 0) {
    matrix(stats::rnorm(n_id * D * K_idm), nrow = n_id * D, ncol = K_idm)
  } else {
    matrix(0.0, nrow = n_id * D, ncol = 0)
  }

  re_idm <- if (K_idm > 0) {
    aperm(array(re_idm_flat, dim = c(D, n_id, K_idm)), c(2, 1, 3))
  } else {
    array(0.0, dim = c(n_id, D, 0))
  }

  M_cov <- if (K_idm > 0) {
    if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) K_idm else K_idm * (K_idm + 1L) / 2L
  } else {
    0L
  }

  idx_row_cov <- integer(M_cov)
  idx_col_cov <- integer(M_cov)
  if (M_cov > 0) {
    if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) {
      idx_row_cov <- seq_len(M_cov)
      idx_col_cov <- seq_len(M_cov)
    } else {
      pos <- 1L
      for (r in seq_len(K_idm)) {
        for (c in seq_len(r)) {
          idx_row_cov[pos] <- r
          idx_col_cov[pos] <- c
          pos <- pos + 1L
        }
      }
    }
  }

  .sim_cov_lp_to_matrix <- function(lp_vec, q_idm, idx_row, idx_col, diag_link) {
    mat <- matrix(0.0, nrow = q_idm, ncol = q_idm)
    if (q_idm <= 0 || length(lp_vec) == 0) {
      return(mat)
    }
    for (m in seq_along(lp_vec)) {
      r_ <- idx_row[m]
      c_ <- idx_col[m]
      val <- as.numeric(lp_vec[m])
      if (r_ == c_) {
        val <- if (diag_link == "exp") exp(val) else log1p(exp(val))
      }
      mat[r_, c_] <- val
    }
    mat
  }

  .sim_cov_matrix_to_alpha <- function(mat, idx_row, idx_col, diag_link) {
    if (length(idx_row) == 0) return(numeric(0))
    out <- numeric(length(idx_row))
    for (m in seq_along(idx_row)) {
      r_ <- idx_row[m]
      c_ <- idx_col[m]
      val <- mat[r_, c_]
      if (r_ == c_) {
        if (!is.finite(val) || val <= 0) {
          cli::cli_abort(c(
            x = "Translated baseline covariance produced a non-positive diagonal entry.",
            i = "Check the supplied {.arg re_params$id_marker_cov$latent} values."
          ))
        }
        out[m] <- .sim_inverse_diag_link(val, diag_link)
      } else {
        out[m] <- val
      }
    }
    out
  }

  .sim_align_len <- function(x, n, default) {
    if (n <= 0) return(numeric(0))
    if (is.null(x)) {
      default <- as.numeric(default)
      if (length(default) == 1L) return(rep(default, n))
      if (length(default) == n) return(default)
      cli::cli_abort(c(
        x = "Internal covariance-regression default has incompatible length.",
        i = "Expected length {n}, got {length(default)}."
      ))
    }
    x <- as.numeric(x)
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) != n) {
      cli::cli_abort(c(
        x = "Length mismatch in covariance-regression parameter specification.",
        i = "Expected length {n}, got {length(x)}."
      ))
    }
    x
  }

  .sim_canonicalize_cov_latent <- function(lambda, z) {
    lambda <- as.numeric(lambda)
    if (length(lambda) == 0L) {
      return(list(
        lambda = numeric(0),
        z = matrix(0.0, nrow = nrow(z), ncol = 0L),
        sign = numeric(0)
      ))
    }

    sign_vec <- ifelse(lambda < 0, -1, 1)
    lambda_out <- abs(lambda)
    z_out <- z[, seq_along(lambda), drop = FALSE]
    z_out <- sweep(z_out, 2L, sign_vec, `*`)

    list(lambda = lambda_out, z = z_out, sign = sign_vec)
  }

  .sim_align_cov_beta <- function(beta_cfg, m_cov, k_cov) {
    if (m_cov <= 0) return(matrix(0.0, 0, k_cov))
    if (k_cov <= 0) return(matrix(0.0, m_cov, 0))
    if (is.null(beta_cfg)) {
      return(matrix(stats::rnorm(m_cov * k_cov, mean = 0, sd = 0.12), nrow = m_cov, ncol = k_cov))
    }
    if (is.matrix(beta_cfg)) {
      if (!all(dim(beta_cfg) == c(m_cov, k_cov))) {
        cli::cli_abort(c(
          x = "{.arg re_params$id_marker_cov$beta} matrix has incompatible dimensions.",
          i = "Expected {m_cov}x{k_cov}, got {nrow(beta_cfg)}x{ncol(beta_cfg)}."
        ))
      }
      beta_mat <- matrix(as.numeric(beta_cfg), nrow = m_cov, ncol = k_cov)
      if (any(!is.finite(beta_mat))) {
        cli::cli_abort(c(
          x = "{.arg re_params$id_marker_cov$beta} must contain only finite values.",
          i = "Check the supplied covariance-regression coefficient matrix."
        ))
      }
      return(beta_mat)
    }
    beta_vec <- as.numeric(beta_cfg)
    if (length(beta_vec) == k_cov) {
      beta_mat <- matrix(rep(beta_vec, each = m_cov), nrow = m_cov, ncol = k_cov, byrow = FALSE)
      if (any(!is.finite(beta_mat))) {
        cli::cli_abort(c(
          x = "{.arg re_params$id_marker_cov$beta} must contain only finite values.",
          i = "Check the supplied covariance-regression coefficients."
        ))
      }
      return(beta_mat)
    }
    if (length(beta_vec) == m_cov * k_cov) {
      beta_mat <- matrix(beta_vec, nrow = m_cov, ncol = k_cov, byrow = TRUE)
      if (any(!is.finite(beta_mat))) {
        cli::cli_abort(c(
          x = "{.arg re_params$id_marker_cov$beta} must contain only finite values.",
          i = "Check the supplied covariance-regression coefficients."
        ))
      }
      return(beta_mat)
    }
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$beta} has incompatible length.",
      i = "Provide length {k_cov}, or {m_cov * k_cov}, or an explicit {m_cov}x{k_cov} matrix."
    ))
  }

  diag_link_cov <- tolower(as.character(re_cov_cfg$diag_link %||% "softplus")[1])
  if (!diag_link_cov %in% c("softplus", "exp")) {
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$diag_link} must be 'softplus' or 'exp'.",
      i = "Use 'softplus' (default) or 'exp'."
    ))
  }

  default_alpha_cov <- if (M_cov > 0L) {
    vapply(seq_len(M_cov), function(m) {
      r_ <- idx_row_cov[m]
      c_ <- idx_col_cov[m]
      if (r_ == c_) {
        diag_target <- stats::runif(1, min = 0.45, max = 1.15)
        if (diag_link_cov == "exp") log(diag_target) else log(expm1(diag_target))
      } else {
        stats::rnorm(1, mean = 0, sd = 0.2)
      }
    }, numeric(1))
  } else {
    numeric(0)
  }
  default_lambda_cov <- if (M_cov > 0L) {
    .sim_random_signed(M_cov, min_abs = 0.15, max_abs = 0.75)
  } else {
    numeric(0)
  }
  alpha_cov <- .sim_align_len(re_cov_cfg$alpha, M_cov, default = default_alpha_cov)
  lambda_cov <- .sim_align_len(re_cov_cfg$lambda, M_cov, default = default_lambda_cov)
  if (any(!is.finite(alpha_cov))) {
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$alpha} must contain only finite values.",
      i = "Check the supplied covariance-regression intercepts."
    ))
  }
  if (any(!is.finite(lambda_cov))) {
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$lambda} must contain only finite values.",
      i = "Check the supplied covariance-regression loadings."
    ))
  }
  if (!is.null(re_cov_cfg$sd_u) || !is.null(re_cov_cfg$tau_u)) {
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$sd_u} is no longer supported.",
      i = "Covariance regression now uses {.arg lambda} on an iid standard-normal latent directly.",
      i = "Remove {.arg sd_u} or {.arg tau_u} and rescale {.arg lambda} instead."
    ))
  }
  if (!is.null(re_idm_legacy) && M_cov > 0 && K_idm > 0) {
    base_li <- .sim_cov_lp_to_matrix(alpha_cov, K_idm, idx_row_cov, idx_col_cov, diag_link_cov)
    alpha_cov <- .sim_cov_matrix_to_alpha(base_li %*% latent_compat_factor, idx_row_cov, idx_col_cov, diag_link_cov)
  }
  beta_cov <- .sim_align_cov_beta(re_cov_cfg$beta, M_cov, K_cov)

  L_i <- array(0.0, dim = c(n_id, K_idm, K_idm))
  z_cov <- matrix(0.0, nrow = n_id, ncol = max(1L, M_cov))
  lambda_cov_sign <- rep(1, M_cov)
  if (M_cov > 0 && K_idm > 0) {
    z_cov_raw <- matrix(stats::rnorm(n_id * M_cov), nrow = n_id, ncol = M_cov)
    cov_latent <- .sim_canonicalize_cov_latent(lambda_cov, z_cov_raw)
    lambda_cov <- cov_latent$lambda
    z_cov <- cov_latent$z
    lambda_cov_sign <- cov_latent$sign

    lp_cov <- matrix(alpha_cov, nrow = n_id, ncol = M_cov, byrow = TRUE)
    if (K_cov > 0) {
      lp_cov <- lp_cov + Xcov %*% t(beta_cov)
    }
    lp_cov <- lp_cov + sweep(z_cov, 2L, lambda_cov, `*`)

    for (m in seq_len(M_cov)) {
      r_ <- idx_row_cov[m]
      c_ <- idx_col_cov[m]
      vals <- lp_cov[, m]
      if (r_ == c_) {
        vals <- if (diag_link_cov == "exp") exp(vals) else log1p(exp(vals))
      }
      L_i[cbind(seq_len(n_id), r_, c_)] <- vals
    }
  }

  re_idm_scaled <- array(0.0, dim = c(n_id, D, K_idm))
  if (K_idm > 0) {
    for (i in seq_len(n_id)) {
      Li <- matrix(L_i[i, , ], K_idm, K_idm)
      re_idm_scaled[i, , ] <- re_idm[i, , ] %*% t(Li)
    }
  }

  # ---- Association terms driving survival (syntax aligned with fit)
  if (!is.null(assoc_from_formula)) assoc <- assoc_from_formula
  assoc <- unique(assoc)
  if (length(assoc) == 0) assoc <- "cv_total"
  assoc <- .validate_assoc_channels(assoc, context = "simulate_joinme()")

  M_corr <- .assoc_cov_feature_count(K_idm, include_diag = FALSE)
  M_vcov <- if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) {
    as.integer(K_idm)
  } else {
    .assoc_cov_feature_count(K_idm, include_diag = TRUE)
  }

  assoc_coef_scalar <- c(
    cv_total = 0,
    cv_mean = 0,
    cv_marker = 0,
    cs_total = 0,
    cs_mean = 0,
    cs_marker = 0
  )
  assoc_coef_corr <- rep(0.0, M_corr)
  assoc_coef_vcov <- rep(0.0, M_vcov)

  .sim_extract_named_assoc_vector <- function(x, prefix) {
    if (length(x) == 0) return(numeric(0))
    nms <- names(x)
    if (is.null(nms)) return(numeric(0))
    idx <- grepl(paste0("^", prefix, "($|\\[[0-9]+\\]$|[._]?[0-9]+$)"), nms)
    if (!any(idx)) return(numeric(0))
    vals <- as.numeric(x[idx])
    nms_sel <- nms[idx]
    ord_key <- rep(NA_integer_, length(nms_sel))
    ord_key[nms_sel == prefix] <- 1L
    idx_num <- nms_sel != prefix
    if (any(idx_num)) {
      ord_key[idx_num] <- suppressWarnings(as.integer(gsub("[^0-9]", "", nms_sel[idx_num])))
      ord_key[is.na(ord_key)] <- seq_len(sum(is.na(ord_key))) + 1L
    }
    vals[order(ord_key)]
  }

  .sim_fill_assoc_coefs <- function(assoc, assoc_coefs, M_corr, M_vcov) {
    scalar <- assoc_coef_scalar
    corr <- rep(0.0, M_corr)
    vcov <- rep(0.0, M_vcov)

    .assign_assoc_vector <- function(target, supplied, expected, channel, source = NULL) {
      vals <- as.numeric(supplied)
      if (expected <= 0L || length(vals) == 0L) {
        return(target)
      }
      take <- min(expected, length(vals))
      if (take > 0L) {
        target[seq_len(take)] <- vals[seq_len(take)]
      }
      if (length(vals) > expected) {
        cli::cli_warn(c(
          x = "Ignoring extra {.arg {channel}} association coefficients in {.fn simulate_joinme}.",
          i = "Model defines {expected} {.val {channel}} component{?s}, but {length(vals)} value{?s} were supplied{if (!is.null(source)) paste0(' via ', source) else ''}."
        ))
      }
      target
    }

    if (is.list(assoc_coefs) && !is.null(assoc_coefs$corr) && M_corr > 0) {
      corr <- .assign_assoc_vector(corr, assoc_coefs$corr, M_corr, "corr", source = "assoc_coefs$corr")
    }
    if (is.list(assoc_coefs) && !is.null(assoc_coefs$vcov) && M_vcov > 0) {
      vcov <- .assign_assoc_vector(vcov, assoc_coefs$vcov, M_vcov, "vcov", source = "assoc_coefs$vcov")
    }

    if (!is.null(names(assoc_coefs))) {
      named_terms <- intersect(names(scalar), names(assoc_coefs))
      named_terms <- intersect(named_terms, assoc)
      if (length(named_terms) > 0) {
        scalar[named_terms] <- as.numeric(assoc_coefs[named_terms])
      }

      if (M_corr > 0 && !is.list(assoc_coefs)) {
        vc_named <- .sim_extract_named_assoc_vector(assoc_coefs, "corr")
        if (length(vc_named) > 0) {
          corr <- .assign_assoc_vector(corr, vc_named, M_corr, "corr", source = "named assoc_coefs")
        }
      }
      if (M_vcov > 0 && !is.list(assoc_coefs)) {
        vc_named <- .sim_extract_named_assoc_vector(assoc_coefs, "vcov")
        if (length(vc_named) > 0) {
          vcov <- .assign_assoc_vector(vcov, vc_named, M_vcov, "vcov", source = "named assoc_coefs")
        }
      }
    } else if (length(assoc_coefs) > 0) {
      vals <- as.numeric(assoc_coefs)
      cursor <- 1L
      for (term in assoc) {
        if (cursor > length(vals)) break
        if (identical(term, "corr")) {
          if (M_corr > 0) {
            remaining <- vals[cursor:length(vals)]
            take <- min(M_corr, length(remaining))
            if (take > 0) {
              corr <- .assign_assoc_vector(corr, remaining, M_corr, "corr", source = "positional assoc_coefs")
              cursor <- cursor + take
            }
          }
        } else if (identical(term, "vcov")) {
          if (M_vcov > 0) {
            remaining <- vals[cursor:length(vals)]
            take <- min(M_vcov, length(remaining))
            if (take > 0) {
              vcov <- .assign_assoc_vector(vcov, remaining, M_vcov, "vcov", source = "positional assoc_coefs")
              cursor <- cursor + take
            }
          }
        } else if (term %in% names(scalar)) {
          scalar[[term]] <- vals[cursor]
          cursor <- cursor + 1L
        }
      }
    }

    list(scalar = scalar, corr = corr, vcov = vcov)
  }

  assoc_coef_parts <- if (assoc_coefs_missing || is.null(assoc_coefs) || length(assoc_coefs) == 0) {
    .sim_random_assoc_parts(assoc, M_corr, M_vcov)
  } else {
    .sim_fill_assoc_coefs(assoc, assoc_coefs, M_corr, M_vcov)
  }
  assoc_coef_scalar <- assoc_coef_parts$scalar
  assoc_coef_corr <- assoc_coef_parts$corr
  assoc_coef_vcov <- assoc_coef_parts$vcov

  weighted_terms <- intersect(c("cv_total", "cs_total", "cv_marker", "cs_marker"), assoc)
  if (length(weighted_terms) > 0) {
    assoc_coef_scalar[weighted_terms] <- abs(assoc_coef_scalar[weighted_terms])
  }

  # Use signed marker weights directly (no logistic bounding), aligned with Stan.

  assoc_coef_vec <- assoc_coef_scalar[intersect(names(assoc_coef_scalar), assoc)]
  if ("corr" %in% assoc && M_corr > 0) {
    vc_names <- paste0("corr[", seq_len(M_corr), "]")
    assoc_coef_vec <- c(assoc_coef_vec, stats::setNames(assoc_coef_corr, vc_names))
  }
  if ("vcov" %in% assoc && M_vcov > 0) {
    vc_names <- paste0("vcov[", seq_len(M_vcov), "]")
    assoc_coef_vec <- c(assoc_coef_vec, stats::setNames(assoc_coef_vcov, vc_names))
  }

  transforms <- .sim_normalize_transform_list(transforms)
  .sim_make_component_transform_set <- function(spec, term_name, n_components) {
    # Expand a single user-facing covariance-style transform spec into a list of
    # per-component evaluators. Each component currently shares the same runtime
    # transform definition in simulation, which mirrors the common-input user
    # layout used when building standata for Stan.
    n_components <- as.integer(n_components %||% 0L)
    if (n_components <= 0L) {
      return(list())
    }
    rep(list(.sim_make_assoc_transform(spec, term_name)), n_components)
  }
  tf_funs <- list(
    cv_total = .sim_make_assoc_transform(transforms$cv_total, "cv_total"),
    cv_mean = .sim_make_assoc_transform(transforms$cv_mean, "cv_mean"),
    cv_marker = .sim_make_assoc_transform(transforms$cv_marker, "cv_marker"),
    cs_total = .sim_make_assoc_transform(transforms$cs_total, "cs_total"),
    cs_mean = .sim_make_assoc_transform(transforms$cs_mean, "cs_mean"),
    cs_marker = .sim_make_assoc_transform(transforms$cs_marker, "cs_marker"),
    corr = .sim_make_component_transform_set(transforms$corr, "corr", M_corr),
    vcov = .sim_make_component_transform_set(transforms$vcov, "vcov", M_vcov)
  )

  has_tf_cv_marker <- !is.null(transforms$cv_marker)
  has_tf_cv_total <- !is.null(transforms$cv_total)

  # ---- Latent trajectory evaluators at arbitrary (id, marker, time)
 
  assoc_row_template <- vector("list", n_id)
  marker_factor <- factor(marker_levels, levels = marker_levels)
  for (i in seq_len(n_id)) {
    row_df_i <- dataEvent[rep(i, D), , drop = FALSE]
    row_df_i[[marker_var]] <- marker_factor
    assoc_row_template[[i]] <- row_df_i
  }

  eta_components_all_markers <- function(i, t) {
    row_df <- assoc_row_template[[i]]
    row_df[[time_var]] <- t

    x_fix <- .mm(fixed_rhs, row_df)
    z_id <- .sim_rhs_matrix(id_rhs_list, row_df)
    z_mk <- .sim_rhs_matrix(mk_rhs_list, row_df)
    z_idm <- .sim_rhs_matrix(idm_rhs_list, row_df)

    fixed_part <- if (ncol(x_fix) > 0) as.numeric(x_fix %*% beta_long) else rep(0, D)
    id_part <- if (ncol(z_id) > 0) as.numeric(z_id %*% re_id[i, ]) else rep(0, D)
    mk_part <- if (ncol(z_mk) > 0) rowSums(z_mk * re_marker) else rep(0, D)
    idm_part <- if (ncol(z_idm) > 0) rowSums(z_idm * re_idm_scaled[i, , ]) else rep(0, D)

    mu_mean <- fixed_part + id_part
    mu_marker <- mk_part + idm_part
    list(
      mu_mean = mu_mean,
      mu_marker = mu_marker,
      mu_total = mu_mean + mu_marker
    )
  }

  .sim_corr_features <- function(i) {
    # Build subject-level raw correlation features in the same lower-triangular
    # ordering used by Stan: (2,1), (3,1), (3,2), ... .
    #
    # Important: this returns raw correlation features only.
    # Association weighting (assoc_coef_corr) is applied later *after* the
    # optional transform, so simulation matches Stan semantics:
    #   sum_j a_corr[j] * transform(corr_raw[j])
    if (K_idm < 2) return(numeric(0))
    Li <- matrix(L_i[i, , ], nrow = K_idm, ncol = K_idm)
    if (!all(is.finite(Li))) return(rep(0.0, M_corr))

    sigma_i <- Li %*% t(Li)
    corr_re <- matrix(0.0, nrow = K_idm, ncol = K_idm)
    for (r in seq_len(K_idm)) {
      var_r <- sigma_i[r, r]
      if (!is.finite(var_r) || var_r <= 0) next
      corr_re[r, r] <- 1.0
      if (r > 1L) {
        for (c in seq_len(r - 1L)) {
          var_c <- sigma_i[c, c]
          if (!is.finite(var_c) || var_c <= 0) next
          denom <- sqrt(var_r * var_c)
          if (!is.finite(denom) || denom <= 0) next
          val <- sigma_i[r, c] / denom
          if (!is.finite(val)) next
          val <- max(-0.999999, min(0.999999, val))
          corr_re[r, c] <- val
          corr_re[c, r] <- val
        }
      }
    }

    out <- numeric(M_corr)
    m <- 1L
    for (r in 2:K_idm) {
      for (c in 1:(r - 1L)) {
        out[m] <- max(-0.999999, min(0.999999, corr_re[r, c]))
        m <- m + 1L
      }
    }
    out
  }

  .sim_vcov_features <- function(i) {
    if (K_idm < 1) return(numeric(0))
    Li <- matrix(L_i[i, , ], nrow = K_idm, ncol = K_idm)
    if (!all(is.finite(Li))) {
      return(rep(0.0, M_vcov))
    }
    .assoc_vcov_features_from_chol(
      Li,
      diagonal_only = as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L
    )
  }

  corr_by_id <- lapply(seq_len(n_id), .sim_corr_features)
  vcov_by_id <- lapply(seq_len(n_id), .sim_vcov_features)

  assoc_components <- function(i, t) {
    # Step A: build marker-resolved CV components at t and t + eps_cs.
    now <- eta_components_all_markers(i, t)
    eps <- eta_components_all_markers(i, t + eps_cs)

    w_mean <- function(x) sum(marker_weights_eff * x) / D

    # Step B: aggregate raw CV summaries (mean, marker, total).
    cv_mean_raw <- mean(now$mu_mean)
    cv_marker_raw <- w_mean(now$mu_marker)
    cv_total_raw <- w_mean(now$mu_total)

    cv_mean <- tf_funs$cv_mean(cv_mean_raw)
    cv_marker <- if (has_tf_cv_marker) w_mean(tf_funs$cv_marker(now$mu_marker)) else cv_marker_raw
    cv_total <- if (has_tf_cv_total) w_mean(tf_funs$cv_total(now$mu_total)) else cv_total_raw

    # Step C: compute mean slope via finite difference.
    cs_mean_raw <- (mean(eps$mu_mean) - cv_mean_raw) / eps_cs

    # Step D: marker-level slopes for marker and total components.
    # We keep these at marker resolution so CS transforms are applied per marker
    # before weighted averaging, mirroring Stan semantics.
    cs_marker_raw_by_marker <- (eps$mu_marker - now$mu_marker) / eps_cs
    cs_total_raw_by_marker <- (eps$mu_total - now$mu_total) / eps_cs

    # Step E: CS aggregation semantics (aligned with CV aggregation semantics):
    # transform first at marker level, then weighted-average across markers.
    cs_marker <- w_mean(tf_funs$cs_marker(cs_marker_raw_by_marker))
    cs_total <- w_mean(tf_funs$cs_total(cs_total_raw_by_marker))

    # Step F: keep weighted raw summaries for debugging/inspection helpers.
    cs_marker_raw <- w_mean(cs_marker_raw_by_marker)
    cs_total_raw <- w_mean(cs_total_raw_by_marker)

    corr_vals <- corr_by_id[[i]]
    # CORR semantics are transform-first, then weight:
    #   corr_assoc = sum_j assoc_coef_corr[j] * tf_corr(corr_vals[j])
    # This mirrors Stan and intentionally avoids tf_corr(sum_j a_j * corr_j).
    corr_vals_tf <- if (length(corr_vals) > 0) {
      vapply(seq_along(corr_vals), function(m) {
        as.numeric(tf_funs$corr[[m]](corr_vals[m]))[1] - as.numeric(tf_funs$corr[[m]](0))[1]
      }, numeric(1))
    } else {
      numeric(0)
    }
    corr_assoc <- if (length(corr_vals_tf) > 0) sum(assoc_coef_corr * corr_vals_tf) else 0
    vcov_vals <- vcov_by_id[[i]]
    # VCOV follows the same calculation layout as CORR:
    #   Step 1: extract raw time-constant Cholesky features for subject i.
    #   Step 2: transform each component separately.
    #   Step 3: apply the corresponding association coefficient after the
    #           transform, never before it.
    vcov_vals_tf <- if (length(vcov_vals) > 0) {
      vapply(seq_along(vcov_vals), function(m) {
        as.numeric(tf_funs$vcov[[m]](vcov_vals[m]))[1] - as.numeric(tf_funs$vcov[[m]](0))[1]
      }, numeric(1))
    } else {
      numeric(0)
    }
    vcov_assoc <- if (length(vcov_vals_tf) > 0) sum(assoc_coef_vcov * vcov_vals_tf) else 0

    list(
      cv_total = cv_total,
      cv_mean = cv_mean,
      cv_marker = cv_marker,
      cs_total = cs_total,
      cs_mean = tf_funs$cs_mean(cs_mean_raw),
      cs_marker = cs_marker,
      corr = corr_assoc,
      vcov = vcov_assoc,
      raw = list(
        cv_mean = cv_mean_raw,
        cs_total = cs_total_raw,
        cs_marker = cs_marker_raw,
        cs_mean = cs_mean_raw,
        corr = corr_assoc,
        corr_vals = corr_vals,
        corr_vals_tf = corr_vals_tf,
        vcov = vcov_assoc,
        vcov_vals = vcov_vals,
        vcov_vals_tf = vcov_vals_tf,
        marker_values = now
      )
    )
  }

  h0_fn <- if (is.function(h0)) {
    h0
  } else {
    .sim_make_baseline_hazard(
      spec = baseline_hazard,
      formula_basehaz = formulaBasehaz,
      beta_basehaz = beta_basehaz
    )
  }

  eta_event_i <- as.numeric(W_event %*% beta_event)

  .sim_assoc_lp_from_components <- function(comp) {
    assoc_coef_scalar[["cv_total"]] * comp$cv_total +
      assoc_coef_scalar[["cv_mean"]] * comp$cv_mean +
      assoc_coef_scalar[["cv_marker"]] * comp$cv_marker +
      assoc_coef_scalar[["cs_total"]] * comp$cs_total +
      assoc_coef_scalar[["cs_mean"]] * comp$cs_mean +
      assoc_coef_scalar[["cs_marker"]] * comp$cs_marker +
      comp$corr +
      comp$vcov
  }

  .sim_time_key <- function(t) sprintf("%.12f", as.numeric(t))
  assoc_lp_cache <- replicate(n_id, new.env(parent = emptyenv(), hash = TRUE), simplify = FALSE)
  cumhaz_cache <- replicate(n_id, new.env(parent = emptyenv(), hash = TRUE), simplify = FALSE)

  .assoc_lp_cached <- function(i, t) {
    key <- .sim_time_key(t)
    cache_i <- assoc_lp_cache[[i]]
    if (exists(key, envir = cache_i, inherits = FALSE)) {
      return(get(key, envir = cache_i, inherits = FALSE))
    }
    lp_val <- .sim_assoc_lp_from_components(assoc_components(i, t))
    assign(key, lp_val, envir = cache_i)
    lp_val
  }

  # ---- Resolve Gauss-Kronrod nodes/weights for survival integration
  gk_spec <- gk_quadrature(nodes = quadrature_nodes)

  .sim_composite_gk <- function(f, upper, panels) {
    # Composite Gauss-Kronrod with nodes defined on [0,1] per panel.
    if (upper <= 0) return(0)
    edges <- seq(0, upper, length.out = panels + 1L)
    total <- 0
    for (p in seq_len(panels)) {
      a <- edges[p]
      b <- edges[p + 1L]
      nodes <- a + (b - a) * gk_spec$nodes
      vals <- f(nodes)
      total <- total + (b - a) * sum(gk_spec$weights * vals)
    }
    total
  }

  hazard_i <- function(i, t) {
    t <- as.numeric(t)
    if (length(t) > 1L) {
      return(vapply(t, function(tt) h0_fn(tt) * exp(eta_event_i[i] + .assoc_lp_cached(i, tt)), numeric(1)))
    }
    h0_fn(t) * exp(eta_event_i[i] + .assoc_lp_cached(i, t))
  }

  cumhaz_i <- function(i, t) {
    t <- as.numeric(t)
    if (length(t) > 1L) return(vapply(t, function(tt) cumhaz_i(i, tt), numeric(1)))
    if (t <= 0) return(0)

    key <- .sim_time_key(t)
    cache_i <- cumhaz_cache[[i]]
    if (exists(key, envir = cache_i, inherits = FALSE)) {
      return(get(key, envir = cache_i, inherits = FALSE))
    }

    method <- tolower(as.character(integration_control$method %||% "integrate")[1])
    rel_tol <- as.numeric(integration_control$rel.tol %||% 1e-6)
    max_panels <- max(1L, as.integer(integration_control$subdivisions %||% 512L))
    max_refine <- max(1L, as.integer(integration_control$max_refine %||% 8L))

    value <- if (identical(method, "integrate")) {
      integrate_ctrl <- integration_control
      integrate_ctrl$method <- NULL
      integrate_ctrl$max_refine <- NULL
      out <- do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integrate_ctrl))
      as.numeric(out$value)
    } else {
      panels <- 1L
      prev <- NA_real_
      curr <- NA_real_
      for (iter in seq_len(max_refine)) {
        curr <- .sim_composite_gk(function(u) hazard_i(i, u), upper = t, panels = panels)
        if (is.finite(prev) && abs(curr - prev) <= rel_tol * max(1, abs(prev))) break
        prev <- curr
        if (panels >= max_panels) break
        panels <- min(max_panels, panels * 2L)
      }
      as.numeric(curr)
    }

    assign(key, value, envir = cache_i)
    value
  }

  draw_event_time <- function(i) {
    U <- stats::runif(1)
    target_H <- -log(U)
    f_root <- function(t) cumhaz_i(i, t) - target_H
    br <- .find_bracket(
      f_root,
      lower = 0,
      upper = root_control$t_init,
      upper_max = root_control$t_max,
      expand = root_control$expand,
      max_expand = root_control$max_expand
    )
    if (is.null(br)) {
      return(list(time = time_cens, event = 0L, bracketing_failed = TRUE))
    }
    T_star <- stats::uniroot(f_root, lower = br$lower, upper = br$upper)$root
    if (T_star > time_cens) {
      list(time = time_cens, event = 0L, bracketing_failed = FALSE)
    } else {
      list(time = T_star, event = 1L, bracketing_failed = FALSE)
    }
  }

  # ---- Draw event/censoring times
  event_draws <- .sim_parallel_lapply(seq_len(n_id), draw_event_time)

  .is_valid_event_draw <- function(x) {
    is.list(x) &&
      !is.null(x$time) && length(x$time) == 1L && is.finite(as.numeric(x$time)) &&
      !is.null(x$event) && length(x$event) == 1L && is.finite(as.numeric(x$event))
  }

  bad_idx <- which(!vapply(event_draws, .is_valid_event_draw, logical(1)))
  if (length(bad_idx) > 0L) {
    cli::cli_warn(c(
      x = "{length(bad_idx)} parallel event-draw job(s) returned malformed output.",
      i = "Recomputing those jobs deterministically in-process."
    ))

    for (idx in bad_idx) {
      set.seed(seed_base + as.integer(idx))
      event_draws[[idx]] <- draw_event_time(idx)
    }

    still_bad <- which(!vapply(event_draws, .is_valid_event_draw, logical(1)))
    if (length(still_bad) > 0L) {
      cli::cli_abort(c(
        x = "Failed to generate valid event draws for indices: {paste(still_bad, collapse = ', ')}.",
        i = "Check custom hazard/transforms for numerical failures."
      ))
    }
  }

  dataEvent[[event_time_var]] <- vapply(event_draws, function(x) as.numeric(x$time), numeric(1))
  dataEvent[[event_var]] <- vapply(event_draws, function(x) as.integer(x$event), integer(1))

  # ---- Build longitudinal observation schedule conditional on event times
  n_obs_target <- max(2L, as.integer(n_obs_per_marker_per_id))
  cov_names <- setdiff(colnames(dataEvent), c(id_var, event_time_var, event_var))

  .sim_obs_rows_for_id <- function(i) {
    obs_upper <- max(1e-8, min(dataEvent[[event_time_var]][i], time_cens))
    rows_i <- vector("list", D)
    for (d in seq_len(D)) {
      if (!is.null(times_obs) && length(times_obs) > 0) {
        candidate_times <- times_obs[times_obs <= obs_upper]
        if (length(candidate_times) == 0) {
          t_obs <- sort(stats::runif(n_obs_target, 0, obs_upper))
        } else {
          t_obs <- sort(sample(candidate_times, size = n_obs_target, replace = length(candidate_times) < n_obs_target))
        }
      } else {
        t_obs <- sort(stats::runif(n_obs_target, 0, obs_upper))
      }

      row_df <- data.frame(
        id = rep(dataEvent[[id_var]][i], length(t_obs)),
        marker = factor(rep(marker_levels[d], length(t_obs)), levels = marker_levels),
        time = t_obs,
        stringsAsFactors = FALSE
      )
      names(row_df)[names(row_df) == "id"] <- id_var
      names(row_df)[names(row_df) == "time"] <- time_var
      names(row_df)[names(row_df) == "marker"] <- marker_var

      for (cov_nm in cov_names) {
        row_df[[cov_nm]] <- dataEvent[[cov_nm]][i]
      }
      rows_i[[d]] <- row_df
    }
    do.call(rbind, rows_i)
  }

  obs_rows <- .sim_parallel_lapply(seq_len(n_id), .sim_obs_rows_for_id)
  dataLong <- do.call(rbind, obs_rows)
  rownames(dataLong) <- NULL

  # ---- Mean structure from model matrices and sampled random effects
  X_long <- .mm(fixed_rhs, dataLong)
  Z_id_long <- .sim_rhs_matrix(id_rhs_list, dataLong)
  Z_mk_long <- .sim_rhs_matrix(mk_rhs_list, dataLong)
  Z_idm_long <- .sim_rhs_matrix(idm_rhs_list, dataLong)

  id_index <- match(as.character(dataLong[[id_var]]), as.character(dataEvent[[id_var]]))
  marker_index <- match(as.character(dataLong[[marker_var]]), marker_levels)

  mu_long <- as.numeric(X_long %*% beta_long)
  if (K_id > 0) {
    mu_long <- mu_long + rowSums(Z_id_long * re_id[id_index, , drop = FALSE])
  }
  if (K_mk > 0) {
    mu_long <- mu_long + rowSums(Z_mk_long * re_marker[marker_index, , drop = FALSE])
  }
  if (K_idm > 0) {
    re_idm_flat_scaled <- matrix(aperm(re_idm_scaled, c(2, 1, 3)), nrow = n_id * D, ncol = K_idm)
    idx_flat <- (id_index - 1L) * D + marker_index
    idm_effect <- re_idm_flat_scaled[idx_flat, , drop = FALSE]
    mu_long <- mu_long + rowSums(Z_idm_long * idm_effect)
  }

  # ---- Step: apply inverse-link to the linear predictor per marker
  # This yields response-scale means/probabilities used in sampling.
  mu_linked <- vapply(seq_len(nrow(dataLong)), function(r) {
    .sim_apply_inv_link_bc(mu_long[r], inv_link_specs[[marker_index[r]]])
  }, numeric(1))

  # ---- Distributional parameters (family defaults + formulaDist overrides)
  family_by_row <- family_names[marker_index]
  attr(dataLong, "joinme_family_by_row") <- family_by_row

  sigma_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "sigma", 1.0), numeric(1))
  nu_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "nu", 4.0), numeric(1))
  phi_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi", 2.0), numeric(1))
  alpha_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "alpha", 0.0), numeric(1))
  phi_beta_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi_beta", 10.0), numeric(1))
  tau_sde_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "tau_sde", 0.5), numeric(1))
  trials_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "trials", 10L), numeric(1))

  dist_formulas <- .normalize_formula_dist(formulaDist)
  family_names_present <- vapply(sort(unique(family_codes)), .family_code_to_name, character(1))
  .validate_dist_formula_scopes(dist_formulas, family_names_present)
  dist_re_effective <- list()
  dist_coef_effective <- list()
  dist_eta_effective <- list()
  dist_design_cols <- list()
  for (param_name in names(dist_formulas)) {
    X_param <- .build_dist_matrix(dist_formulas[[param_name]], dataLong, family_by_row = family_by_row)$X
    coef_spec <- .sim_dist_coef_spec(param_name, dist_coefs, dist_formulas[[param_name]])
    beta_param <- .sim_align_dist_coef(
      col_names = colnames(X_param),
      user_coef = coef_spec,
      dist_scope = dist_formulas[[param_name]],
      sd_default = 0.15,
      intercept_default = 0.0
    )
    eta_param <- as.numeric(X_param %*% beta_param)

    re_terms_param <- .build_dist_re_terms(dist_formulas[[param_name]], dataLong)
    if (!is.null(re_terms_param$n_re) && re_terms_param$n_re > 0) {
      cfg_dist <- (re_params$dist %||% list())[[param_name]] %||% list()
      dist_re_effective_terms <- vector("list", re_terms_param$n_re)
      for (j in seq_len(re_terms_param$n_re)) {
        cfg_j <- if (!is.null(cfg_dist$terms) && length(cfg_dist$terms) >= j) cfg_dist$terms[[j]] else cfg_dist
        resolved_j <- .sim_resolve_re_block_cfg(
          cfg = cfg_j,
          K = as.integer(re_terms_param$K[j]),
          label = paste0("dist_", param_name, "_re", j)
        )
        dist_re_effective_terms[[j]] <- resolved_j
        b_j <- .sim_draw_re_block(
          n_group = as.integer(re_terms_param$G[j]),
          K = as.integer(re_terms_param$K[j]),
          cfg = cfg_j,
          label = paste0("dist_", param_name, "_re", j),
          resolved = resolved_j
        )
        eta_param <- eta_param + rowSums(re_terms_param$Z[[j]] * b_j[re_terms_param$J[[j]], , drop = FALSE])
      }
      dist_re_effective[[param_name]] <- list(terms = dist_re_effective_terms)
    }

    dist_coef_effective[[param_name]] <- beta_param
    dist_eta_effective[[param_name]] <- eta_param
    dist_design_cols[[param_name]] <- colnames(X_param) %||% character(0)

    if (param_name == "sigma") sigma_vec <- exp(eta_param)
    if (param_name == "nu") nu_vec <- 2 + exp(eta_param)
    if (param_name == "phi") phi_vec <- exp(eta_param)
    if (param_name == "alpha") alpha_vec <- eta_param
    if (param_name == "phi_beta") phi_beta_vec <- exp(eta_param)
    if (param_name == "tau_sde") tau_sde_vec <- stats::plogis(eta_param)
  }

  # ---- Step: draw outcomes by marker-specific family
  y_out <- numeric(nrow(dataLong))
  for (r in seq_len(nrow(dataLong))) {
    fam_name <- family_by_row[r]
    fam_code <- .parse_family(fam_name)
    eta_r <- mu_long[r]
    mu_r <- mu_linked[r]

    if (fam_code == 11L) {
      cp <- family_params[[fam_name]]$cutpoints %||% c(-1, 1)
      y_out[r] <- .sim_sample_ordinal(eta = eta_r, cutpoints = cp)
    } else {
      y_out[r] <- .sample_from_family(
        n = 1,
        mu = mu_r,
        family = fam_code,
        sigma = sigma_vec[r],
        nu = nu_vec[r],
        phi = phi_vec[r],
        phi_beta = phi_beta_vec[r],
        tau_sde = tau_sde_vec[r],
        trials = as.integer(round(trials_vec[r])),
        skew = alpha_vec[r]
      )
    }
  }
  dataLong[[y_var]] <- y_out

  # ---- Final formatting and metadata
  dataLong[[marker_var]] <- factor(dataLong[[marker_var]], levels = marker_levels)
  ord_long <- order(dataLong[[id_var]], dataLong[[marker_var]], dataLong[[time_var]])
  dataLong <- dataLong[ord_long, , drop = FALSE]
  family_by_row <- family_by_row[ord_long]
  sigma_vec <- sigma_vec[ord_long]
  nu_vec <- nu_vec[ord_long]
  phi_vec <- phi_vec[ord_long]
  alpha_vec <- alpha_vec[ord_long]
  phi_beta_vec <- phi_beta_vec[ord_long]
  tau_sde_vec <- tau_sde_vec[ord_long]
  trials_vec <- trials_vec[ord_long]
  rownames(dataLong) <- NULL

  dist_param_rowwise <- data.frame(
    id = dataLong[[id_var]],
    marker = dataLong[[marker_var]],
    time = dataLong[[time_var]],
    family = family_by_row,
    sigma = sigma_vec,
    nu = nu_vec,
    phi = phi_vec,
    phi_nb = phi_vec,
    alpha = alpha_vec,
    alpha_skew = alpha_vec,
    skew = alpha_vec,
    phi_beta = phi_beta_vec,
    tau_sde = tau_sde_vec,
    trials = as.integer(round(trials_vec)),
    stringsAsFactors = FALSE
  )

  distributional_truth <- list(
    family_by_row = family_by_row,
    rowwise = dist_param_rowwise,
    sigma = sigma_vec,
    nu = nu_vec,
    phi = phi_vec,
    phi_nb = phi_vec,
    alpha = alpha_vec,
    alpha_skew = alpha_vec,
    skew = alpha_vec,
    phi_beta = phi_beta_vec,
    tau_sde = tau_sde_vec,
    trials = as.integer(round(trials_vec)),
    coef = dist_coef_effective,
    eta = dist_eta_effective,
    design_cols = dist_design_cols,
    family_defaults = family_params
  )

  tmax <- max(dataEvent[[event_time_var]])

  beta_eff_in_likelihood <- beta_long
  if (length(idx_time_beta) > 0L) {
    beta_eff_in_likelihood[idx_time_beta] <- beta_eff_in_likelihood[idx_time_beta] * tmax
  }

  tau_u_eff <- re_id_effective$sd
  if (length(idx_time_uid) > 0L) {
    tau_u_eff[idx_time_uid] <- tau_u_eff[idx_time_uid] * tmax
  }

  tau_v_eff <- re_marker_effective$sd
  if (length(idx_time_vmk) > 0L) {
    tau_v_eff[idx_time_vmk] <- tau_v_eff[idx_time_vmk] * tmax
  }

  marker_id_row_scale_eff <- rep(1.0, K_idm)
  if (length(idx_time_widm) > 0L) {
    marker_id_row_scale_eff[idx_time_widm] <- tmax
  }

  family_codes_present <- sort(unique(family_codes))
  family_names_present <- vapply(family_codes_present, .family_code_to_name, character(1))

  .sim_family_param_truth <- function(param_name, fallback, include_families = NULL) {
    fam_keep <- family_codes_present[vapply(
      family_codes_present,
      function(fc) param_name %in% .family_distrib_params(fc),
      logical(1)
    )]
    if (!is.null(include_families)) {
      fam_keep <- fam_keep[vapply(fam_keep, function(fc) .family_code_to_name(fc) %in% include_families, logical(1))]
    }
    if (length(fam_keep) == 0L) return(numeric(0))
    fam_names_keep <- vapply(fam_keep, .family_code_to_name, character(1))
    vals <- vapply(fam_names_keep, function(fam) .sim_get_family_param(fam, param_name, fallback), numeric(1))
    names(vals) <- fam_names_keep
    vals
  }

  sigma_family_truth <- .sim_family_param_truth("sigma", 1.0)
  nu_family_truth <- .sim_family_param_truth("nu", 4.0)
  phi_family_truth <- .sim_family_param_truth("phi", 2.0)
  alpha_family_truth <- .sim_family_param_truth("alpha", 0.0)
  phi_beta_family_truth <- .sim_family_param_truth("phi_beta", 10.0)
  tau_sde_family_truth <- .sim_family_param_truth("tau_sde", 0.5)
  trials_family_truth <- if (any(family_names_present == "binomial")) {
    vals <- c(binomial = .sim_get_family_param("binomial", "trials", 10L))
    as.numeric(setNames(vals, names(vals)))
  } else {
    numeric(0)
  }
  cutpoints_ord_truth <- if (any(family_names_present == "cumulative_logit")) {
    as.numeric(family_params$cumulative_logit$cutpoints %||% c(-1, 1))
  } else {
    numeric(0)
  }

  marker_to_sigma_family <- setNames(match(family_names, names(sigma_family_truth), nomatch = 0L), marker_levels)
  marker_to_nu_family <- setNames(match(family_names, names(nu_family_truth), nomatch = 0L), marker_levels)
  marker_to_phi_family <- setNames(match(family_names, names(phi_family_truth), nomatch = 0L), marker_levels)
  marker_to_alpha_family <- setNames(match(family_names, names(alpha_family_truth), nomatch = 0L), marker_levels)
  marker_to_phi_beta_family <- setNames(match(family_names, names(phi_beta_family_truth), nomatch = 0L), marker_levels)
  marker_to_tau_sde_family <- setNames(match(family_names, names(tau_sde_family_truth), nomatch = 0L), marker_levels)

  stan_fit_truth <- list(
    beta = beta_long,
    beta_eff_in_likelihood = beta_eff_in_likelihood,
    gamma_w = beta_event,
    bs_gamma_c = beta_basehaz,
    tau_u = re_id_effective$sd,
    tau_u_eff = tau_u_eff,
    Lcorr_u = re_id_effective$Lcorr,
    Corr_u = re_id_effective$corr,
    Sigma_u = re_id_effective$cov,
    tau_v = re_marker_effective$sd,
    tau_v_eff = tau_v_eff,
    Lcorr_v = re_marker_effective$Lcorr,
    Corr_v = re_marker_effective$corr,
    Sigma_v = re_marker_effective$cov,
    marker_id_row_scale_eff = marker_id_row_scale_eff,
    alpha_cv_total_eff = unname(assoc_coef_scalar[["cv_total"]] %||% NA_real_),
    alpha_cs_total_eff = unname(assoc_coef_scalar[["cs_total"]] %||% NA_real_),
    alpha_cv_marker_eff = unname(assoc_coef_scalar[["cv_marker"]] %||% NA_real_),
    alpha_cs_marker_eff = unname(assoc_coef_scalar[["cs_marker"]] %||% NA_real_),
    alpha_corr_eff = assoc_coef_corr,
    alpha_vcov_eff = assoc_coef_vcov,
    sigma_family = sigma_family_truth,
    nu_family = nu_family_truth,
    phi_family = phi_family_truth,
    alpha_family = alpha_family_truth,
    phi_beta_family = phi_beta_family_truth,
    tau_sde_family = tau_sde_family_truth,
    trials_family = trials_family_truth,
    cutpoints_ord = cutpoints_ord_truth,
    marker_to_sigma_family = marker_to_sigma_family,
    marker_to_nu_family = marker_to_nu_family,
    marker_to_phi_family = marker_to_phi_family,
    marker_to_alpha_family = marker_to_alpha_family,
    marker_to_phi_beta_family = marker_to_phi_beta_family,
    marker_to_tau_sde_family = marker_to_tau_sde_family
  )

  family_truth <- list(
    by_marker = data.frame(
      marker = marker_levels,
      family = family_names,
      family_code = family_codes,
      marker_to_sigma_family = unname(marker_to_sigma_family),
      marker_to_nu_family = unname(marker_to_nu_family),
      marker_to_phi_family = unname(marker_to_phi_family),
      marker_to_alpha_family = unname(marker_to_alpha_family),
      marker_to_phi_beta_family = unname(marker_to_phi_beta_family),
      marker_to_tau_sde_family = unname(marker_to_tau_sde_family),
      stringsAsFactors = FALSE
    ),
    shared = stan_fit_truth[c(
      "sigma_family",
      "nu_family",
      "phi_family",
      "alpha_family",
      "phi_beta_family",
      "tau_sde_family",
      "trials_family",
      "cutpoints_ord"
    )],
    defaults = family_params
  )

  true_params <- list(
    beta_long = beta_long,
    beta_event = beta_event,
    alpha_cv_total = if ("cv_total" %in% names(assoc_coef_vec)) assoc_coef_vec[["cv_total"]] else NA_real_,
    assoc = assoc,
    assoc_coefs = assoc_coef_vec,
    marker_weights = marker_weights_eff,
    marker_weights_raw = marker_weights_raw,
    marker_weights_eff = marker_weights_eff,
    link_names = link_names,
    quadrature_nodes = as.integer(gk_spec$n_gk),
    gk_nodes = gk_spec$nodes,
    gk_weights = gk_spec$weights,
    gk_rule = gk_spec$rule,
    transforms = transforms,
    baseline_hazard = baseline_hazard,
    formulaVCov = formulaVCov,
    formulaBasehaz = formulaBasehaz,
    beta_basehaz = beta_basehaz,
    family = family_truth,
    distributional_params = distributional_truth,
    stan_fit = stan_fit_truth,
    dist_coefs = dist_coefs,
    dist_re_params = re_params$dist %||% list(),
    dist_re_effective = dist_re_effective %||% list(),
    re_params = re_params,
    re_effective = list(
      id = re_id_effective,
      marker = re_marker_effective,
      id_marker_cov = list(
        latent = re_idm_effective,
        alpha = alpha_cov,
        beta = beta_cov,
        lambda = lambda_cov,
        z = z_cov[, seq_len(M_cov), drop = FALSE],
        lambda_sign = lambda_cov_sign,
        diag_link = diag_link_cov
      )
    ),
    re_draws = list(
      id = re_id,
      marker = re_marker,
      id_marker_cov_latent = re_idm,
      id_marker_cov_scaled = re_idm_scaled
    ),
    L_i = L_i,
    id_marker_cov_effective = list(
      latent = re_idm_effective,
      legacy_latent_translation = list(
        input = re_idm_legacy,
        factor = latent_compat_factor,
        applied_to = if (!is.null(re_idm_legacy)) "alpha" else NULL
      ),
      alpha = alpha_cov,
      beta = beta_cov,
      lambda = lambda_cov,
      z = z_cov[, seq_len(M_cov), drop = FALSE],
      lambda_sign = lambda_cov_sign,
      diag_link = diag_link_cov
    ),
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaDist = dist_formulas
  )

  marker_info <- list(
    names = marker_levels,
    families = families,
    family_codes = family_codes
  )

  helpers <- list(
    cv_total = function(i, t) assoc_components(i, t)$cv_total,
    cv_mean = function(i, t) assoc_components(i, t)$cv_mean,
    cv_marker = function(i, t) assoc_components(i, t)$cv_marker,
    cs_total = function(i, t) assoc_components(i, t)$cs_total,
    cs_mean = function(i, t) assoc_components(i, t)$cs_mean,
    cs_marker = function(i, t) assoc_components(i, t)$cs_marker,
    corr = function(i, t) assoc_components(i, t)$corr,
    vcov = function(i, t) assoc_components(i, t)$vcov,
    assoc_components_raw = function(i, t) assoc_components(i, t)$raw,
    assoc_components = assoc_components,
    baseline_hazard = h0_fn,
    hazard = hazard_i,
    cumhaz = cumhaz_i,
    survival_prob = function(i, t) exp(-cumhaz_i(i, t))
  )

  list(
    dataLong = dataLong,
    dataEvent = dataEvent,
    truth = true_params,
    true_params = true_params,
    marker_info = marker_info,
    helpers = helpers,
    tmax = tmax
  )
}
