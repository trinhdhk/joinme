#' @name JoiNMe_simulation
#' @title Simulation for JoiNMe
#'
#' @importFrom stats rnorm rexp runif plogis rbeta rbinom rnbinom rpois rt integrate uniroot
#'
#' @description
#' Simulator for the JoiNMe joint model aligned with the v11 Stan semantics.
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

#' Softplus function
#' @description Numerically stable implementation of softplus transformation
#' 
#' @details
#' \deqn{\operatorname{softplus}(x) = \log(1+\exp(x))}
#' 
#' @param x Numeric input
#' @export
softplus <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' Draw standardized shrinkage latents
#'
#' @param n Number of draws.
#' @param shrinkage Integer family code: 0 Student-t(6), 1 Laplace, 2 Normal.
#' @return Numeric vector of standardized draws.
#' @keywords internal
#' @noRd
.sim_draw_standard_shrinkage <- function(n, shrinkage) {
  n <- as.integer(n)
  if (n <= 0L) return(numeric(0))
  if (shrinkage == 1L) {
    # Inverse-CDF construction for double_exponential(0, 1).
    return(sample(c(-1, 1), n, replace = TRUE) * stats::rexp(n, rate = 1))
  }
  if (shrinkage == 2L) return(stats::rnorm(n))
  stats::rt(n, df = 6)
}

#' Declare a time-varying simulation covariate generator
#'
#' @description
#' Helper used inside `covariate_formulas` in `simulate_joinme()` to create
#' subject-specific stepwise covariates over time.
#' Elsewhere, this function does not evaluate anything.
#' Only a specification list is returned.
#'
#' @param fun Random generator function (for example `rnorm`, `rbinom`).
#' @param n_step Integer vector of length 1 or 2.
#'   - length 1: each subject receives exactly `n_step` step periods.
#'   - length 2: each subject draws its number of periods uniformly from
#'     `min(n_step):max(n_step)`.
#'   Ignored when `steps` is supplied.
#' @param steps Optional vector of positive integer time points where the
#'   step changes occur exactly. Defaults to `NULL`.
#'   Values greater than subject-specific stop times are ignored.
#' @param ... Additional arguments forwarded to `fun`.
#'
#' @return Specification fed into `simulate_joinme()`.
#' @export
time_varyring <- function(fun, n_step = NULL, steps = NULL, ...) {
  structure(
    list(fun = fun, n_step = n_step, steps = steps, args = list(...)),
    class = "JoiNMe_time_varyring_spec"
  )
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
#' @noRd
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
  seed = .Random.seed[[1]],
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
      L_i[i, q, q] <- softplus(lp)
    }
  }

  # w_idm[i,d] = L_i[i] %*% z_w[i,d]
  w_idm <- array(0.0, dim = c(n_id, D, Q_idm))
  for (i in seq_len(n_id)) {
    Li <- matrix(L_i[i, , ], Q_idm, Q_idm)
    for (d in seq_len(D)) w_idm[i, d, ] <- as.numeric(Li %*% z_w[i, d, ])
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
      idm_part <- if (Q_idm == 1) w_idm[i, d, 1] else w_idm[i, d, 1] + w_idm[i, d, 2] * t
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
        mu_idm <- sum(idm_row(tt) * w_idm[i, d, ])
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
    w_idm = w_idm
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

#' Resolve marker-specific family declarations for simulation
#'
#' @description
#' Converts character and structured family declarations into the
#' marker-aligned codes, inverse-link instructions, and fixed skew-Laplace
#' quantiles used by `simulate_joinme()`. The same family parser and fixed-`tau`
#' normaliser are used by fitting, which prevents simulation and estimation
#' from assigning different meanings to a family declaration.
#'
#' @param families Character vector or list containing one family declaration
#'   per marker.
#' @param D Positive integer number of longitudinal markers.
#'
#' @return A list containing marker-aligned `family_codes`, `link_names`,
#'   `inv_link_specs`, `use_tau_fixed`, and `tau_fixed`.
#' @keywords internal
#' @noRd
.sim_resolve_family_specs <- function(families, D) {
  if (.is_single_family_spec(families)) {
    families <- rep(list(families), D)
  } else if (length(families) == 1L && !is.list(families)) {
    families <- rep(list(families), D)
  } else if (!is.list(families)) {
    families <- as.list(families)
  }

  if (length(families) != D) {
    cli::cli_abort(
      "families must have length {D} (number of markers), got {length(families)}"
    )
  }

  specifications <- lapply(families, .extract_family_and_link)
  family_codes <- vapply(
    specifications,
    function(specification) specification$family_code,
    integer(1)
  )
  fixed_tau <- .normalise_fixed_tau_by_marker(
    use_tau_fixed = as.integer(vapply(
      specifications,
      function(specification) !is.na(specification$tau_fixed),
      logical(1)
    )),
    tau_fixed = vapply(
      specifications,
      function(specification) {
        if (is.na(specification$tau_fixed)) 0.5 else specification$tau_fixed
      },
      numeric(1)
    ),
    family_codes = family_codes
  )
  link_names <- vapply(
    specifications,
    function(specification) specification$link_name,
    character(1)
  )
  link_names[is.na(link_names)] <- "custom"

  list(
    family_codes = as.integer(family_codes),
    link_names = link_names,
    inv_link_specs = lapply(
      specifications,
      function(specification) specification$inv_link_bc
    ),
    use_tau_fixed = fixed_tau$use_tau_fixed,
    tau_fixed = fixed_tau$tau_fixed
  )
}

#' Simulate joint model data
#' 
#' @description
#' Simulator for `JoiNMe` that mirrors the fitting syntax as closely as
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
#'   `covariate_formulas`). In nested marker terms, outer `( ... || marker )`
#'   keeps marker-only and marker-by-id blocks independent, while inner
#'   `( ... || id )` makes the marker-by-id covariance diagonal.
#' @param formulaEvent Event/survival formula (same role as in `joinme()`).
#' @param formulaVCov Optional covariance-regression formula for the id-specific
#'   marker-by-id covariance factor (same role as in `joinme_standata()`).
#'   This formula is evaluated on event-level covariates (one row per subject),
#'   must not include random-effect bars `( ... | ... )`, and must not include
#'   the longitudinal time variable.
#'
#'   Internally, this formula drives subject-specific standard deviations and
#'   Cholesky-correlation-factor rows used to build `L_i = SD_i * K_i`; see
#'   `re_params$id_marker_cov`.
#'   The default `~ 1` is supported and gives an intercept-only covariance regression.
#' @param formulaDist Optional distributional regression formulas (same role as in `joinme()`).
#'   Supported LHS parameters are `sigma`, `nu`, `phi`, `alpha` (aliases:
#'   `alpha_skew`, `skew`), `kappa`, and `tau`.
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
#'   `cbrt`, `softplus`/`log1p_exp`, trigonometric and hyperbolic functions,
#'   the standard normal CDF (`Phi`, `pnorm`), and the standard normal
#'   quantile (`inv_Phi`, `qnorm`, `probit`).
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
#'     Simulation deliberately treats these as fixed interpolation pairs. This
#'     differs from model fitting, where knot ordinates are estimated as an
#'     ordered simplex construction and `y` no longer fixes the curve.
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
#' @param marker_weights Optional base marker weights used for association
#'   aggregation. If `shared_marker_weights = TRUE`, supply one numeric vector to
#'   be shared across all active weighted marker-based association terms. If
#'   `shared_marker_weights = FALSE`, you may instead supply a named list with
#'   entries `cv_total`, `cs_total`, `cv_marker`, and `cs_marker`.
#' @param shared_marker_weights Logical. If `TRUE`, all active weighted
#'   marker-based association terms share one marker-weight structure. If
#'   `FALSE`, each active weighted marker-based association term uses its own
#'   marker-weight structure.
#' @param fixed_marker_weights Logical. This has the same meaning as in
#'   `joinme()`. If `TRUE`, `marker_weights` are the effective weights and no
#'   latent perturbation is drawn. If `FALSE` (the default), `marker_weights`
#'   are base weights and the effective weights are
#'   `marker_weights + z_marker_weights`. When base weights are omitted they
#'   default to zero, exactly as they do in `joinme()` when marker weights are
#'   estimated.
#' @param shrinkage Integer selecting the distribution of each standardized
#'   latent marker-weight perturbation when `fixed_marker_weights = FALSE`:
#'   `0` draws Student-t with 6 degrees of freedom, `1` draws standard Laplace,
#'   and `2` draws standard Normal. This is the same switch used by the Stan
#'   priors. It does not change explicitly supplied fixed weights.
#' @param n_id Number of subjects.
#' @param families Marker-specific family names or `jm_family()` declarations.
#'   Use `jm_family()` entries to supply custom `link`/`inv_link` expressions
#'   or to fix the skew-Laplace quantile for a marker, for example
#'   `jm_family("skew_laplace", tau = 0.8)`.
#'   A named probit link applies `Phi` as its inverse link. In formula
#'   expressions, `Phi`/`pnorm` are the standard normal CDF and
#'   `inv_Phi`/`qnorm`/`probit` are the standard normal quantile. Simulation
#'   evaluates the same bytecode instructions and fixed-quantile selection as
#'   fitting and prediction.
#' @param marker_levels Optional marker names; defaults to `m1`, `m2`, ...
#' @param times_obs Scheduled observation time grid used for every `(id, marker)`
#'   before optional visit-time jitter and post-event censoring are applied.
#' @param obs_time_noise_sd Optional Gaussian noise SD added independently to each
#'   simulated observation time after the visit schedule is chosen; noisy visit
#'   times are clipped to `[0, time_cens]`.
#' @param censor_longitudinal_after_event Logical; if TRUE (default), simulated
#'   longitudinal observations are truncated at the subject event time. If FALSE,
#'   the longitudinal schedule may continue after the event time up to
#'   `time_cens`.
#' @param ... Additional unused compatibility arguments. Legacy observation-count
#'   inputs are ignored; the number of scheduled observations is inferred from
#'   `times_obs`.
#' @param seed RNG seed.
#' @param covariate_formulas Named or LHS formulas used to generate event-level covariates,
#'   e.g. `list(x1 ~ rnorm(n_id), x2 ~ rt(n_id, df = 5))`.
#'   Time-varying step covariates are supported via
#'   `time_varyring(fun, n_step, steps, ...)`, for example
#'   `list(x ~ time_varyring(rnorm, c(3, 6), mean = 0, sd = 1))` or
#'   `list(x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1))`.
#'   In this mode, each id receives stepwise periods over `[0, time_cens]`, and
#'   each period value is sampled from `fun`.
#' @param assoc Association components (same names as fit): `cv_total`, `cv_mean`,
#'   `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.
#'
#'   Meaning of each channel:
#'   - `cv_*`: current-value channels from longitudinal trajectories,
#'   - `cs_*`: current-slope channels from finite differences (`eps_cs`),
#'   - `corr`: off-diagonal entries of the subject-specific Cholesky-correlation
#'     factor `K` from marker-by-id random effects.
#'   - `vcov`: the same off-diagonal `K` entries together with the
#'     subject-specific standard deviations.
#'
#'   You can also provide `formulaAssoc = ~ ...` to select channels; when present,
#'   it overrides `assoc`.
#' @param assoc_coefs Association coefficients for hazard terms.
#'   Non-`corr`/`vcov` terms accept scalar values. The `corr` term accepts a vector of
#'   off-diagonal `K` coefficients ordered as `(2,1), (3,1), (3,2), ...` in
#'   lower-triangular row-major order of the marker-by-id random-effect
#'   covariance dimension. The `vcov` term accepts the same off-diagonal `K`
#'   entries followed by the subject-specific standard deviations,
#'   `(2,1), (3,1), (3,2), ..., sd_1, sd_2, ...`; when `||` is used in the
#'   marker-by-id random-effects block, only the standard deviation entries are used.
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
#'   parameters (`sigma`, `nu`, `phi`, `alpha`, `kappa`, `tau`).
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
#'   - `alpha`: baseline linear predictors for the subject-specific covariance
#'     regression entries. Diagonal positions control standard deviations;
#'     off-diagonal positions control the row-wise correlation-factor regression,
#'     on the tanh scale.
#'   - `beta`: covariate effects from `formulaVCov` design matrix (systematic
#'     subject-to-subject covariance shifts by observed covariates),
#'   - `lambda`: a scalar or length-\(M\) numeric vector of loadings. A scalar
#'     is repeated over all \(M\) covariance coordinates; a vector supplies one
#'     loading \(\lambda_m\) for each packed lower-triangular coordinate.
#'     Each loading multiplies only its matching iid standard-normal subject
#'     perturbation, so the collective operation is
#'     `diag(lambda) %*% z_i`, not a full loading matrix. During fitting every
#'     \(\lambda_m\) is constrained to be non-negative. If a negative value is
#'     supplied for simulation, its sign is folded into the corresponding
#'     latent draw, leaving the generated model unchanged whilst following the
#'     identified fitting convention,
#'   - `diag_link`: link for the subject-specific standard deviations (`"softplus"`
#'     or `"exp"`).
#'
#'   Element-wise covariance-regression form is:
#'   `eta_{im} = alpha_m + x_i^T beta_m + lambda_m z_{im}`, with
#'   `lambda_m >= 0` and `z_{im} ~ Normal(0, 1)`. Diagonal entries apply
#'   `diag_link` to give positive subject-specific standard deviations. Off-diagonal
#'   entries are mapped through `tanh` and then assembled row by row into a valid
#'   subject-specific Cholesky-correlation factor `K_i`. The covariance factor is
#'   reconstructed as `L_i = SD_i * K_i`. Marker-by-id latent seeds are sampled as
#'   iid standard normal values and the final marker-by-id effects are obtained as
#'   `b_id = L_i z_id`.
#'
#'   Therefore, in the ordinary non-mixture model, \(\lambda_m\) is the
#'   conditional standard deviation of the unexplained subject heterogeneity
#'   on covariance-predictor coordinate \(m\), before applying `diag_link` or
#'   `tanh`. It is not itself an entry of `L_i`, a covariance, or a correlation.
#'   With `class_type = "corr"` or `"vcov"`, the selected `z_{im}` has a
#'   class-specific location and scale. Conditional on class \(g\), its
#'   contribution to `eta_{im}` consequently has location
#'   `lambda_m * mix_location[g, m]` and distributional scale
#'   `lambda_m * mix_scale[g, m]`.
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
#' @param family_params Family-specific constants used when the corresponding
#'   parameter has no distributional regression. In particular, specify the
#'   Beta mean/sample-size model with \code{beta = list(kappa = ...)} and the
#'   skew double exponential model with
#'   \code{skew_double_exponential = list(sigma = ..., tau = ...)}.
#'   \code{kappa} must be positive. \code{tau} must lie in \eqn{(0,1)}, with
#'   \eqn{0.5} giving the symmetric double exponential distribution.
#'   A marker-specific `tau` in `jm_family()` takes precedence over this shared
#'   family constant and over a `tau` distributional regression for that
#'   marker, matching the fitted likelihood.
#' @param h0 Optional baseline hazard function `h0(t)`
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
#' @param .mixture_specification Private latent-progress simulation
#'   specification assembled by [simulate_joinme_mix()]. Users should call
#'   [simulate_joinme_mix()] rather than supplying this argument directly.
#' @param left_truncation_max Optional non-negative delayed-entry bound. When > 0,
#'   each subject receives a sampled entry in
#'   `[0, min(left_truncation_max, stop_time))`.
#' @param truncate_longitudinal_before_entry Logical; when TRUE and delayed entry
#'   is active (`left_truncation_max > 0`), longitudinal rows observed before the
#'   sampled entry time are removed.
#'
#' @return A list containing `dataLong`, `dataEvent`, `truth` (also available as
#'   `true_params`), `marker_info`, `helpers`, and `tmax`. Marker-weight truth
#'   distinguishes base, latent, and effective values. Baseline-hazard truth is
#'   stored in `truth$baseline_hazard`; when a log-linear representation exists,
#'   its resolved coefficients are also in `truth$stan_fit$bs_gamma_c`.
#'   Marker-specific fixed skew-Laplace quantiles are recorded in
#'   `truth$stan_fit$use_tau_fixed` and `truth$stan_fit$tau_fixed`; realised
#'   row-specific quantiles are recorded in
#'   `truth$distributional$rowwise$tau`.
#'   Hazard-scale association coefficients are stored as `alpha_cv_total`,
#'   `alpha_cs_total`, `alpha_cv_mean`, `alpha_cs_mean`, `alpha_corr`, and
#'   `alpha_vcov`, matching the fitted posterior output names.
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
#' 
#' sim_split <- simulate_joinme(
#'   n_id = 20,
#'   formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
#'   left_truncation_max = 1.0,
#'   seed = 99
#' )
#'
#' sim_tv <- simulate_joinme(
#'   n_id = 20,
#'   formulaEvent = survival::Surv(time_start, time_stop, event) ~ x + x2,
#'   covariate_formulas = list(
#'     x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1),
#'     x2 ~ rnorm(n_id)
#'   ),
#'   left_truncation_max = 0.5
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
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  ...,
  seed = .Random.seed[[1]],
  covariate_formulas = list(
    x1 ~ rnorm(n_id),
    x2 ~ rnorm(n_id)
  ),
  marker_weights = NULL,
  shared_marker_weights = TRUE,
  fixed_marker_weights = FALSE,
  shrinkage = 0L,
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
    skew_double_exponential = list(sigma = 1.0, tau = 0.5),
    beta = list(kappa = 10),
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
  event_var = "event",
  .mixture_specification = NULL
) {
  # Algorithm overview (statistical simulation workflow):
  # Step 1. Parse formulas and design structures to mirror the fitted-model
  #         likelihood parameterisation exactly.
  # Step 2. Simulate subject-level exogenous covariates, including optional
  #         piecewise-constant time-varying processes.
  # Step 3. Draw fixed/random effects and covariance-regression parameters,
  #         then build latent longitudinal trajectories.
  # Step 4. Construct the event hazard from baseline + covariate + association
  #         channels and sample event times by inverse-transform sampling.
  # Step 5. Generate delayed-entry times (left truncation), construct event
  #         interval rows automatically when required by the data-generating
  #         process, and sample longitudinal observations.
  # Step 6. Draw marker responses from family-specific observation models,
  #         collect truth objects, and return simulation outputs.
  set.seed(seed)
  assoc_coefs_missing <- missing(assoc_coefs)
  if (!is.logical(fixed_marker_weights) || length(fixed_marker_weights) != 1L || is.na(fixed_marker_weights)) {
    cli::cli_abort(c(
      x = "{.arg fixed_marker_weights} must be TRUE/FALSE.",
      i = "Use FALSE to simulate Stan's latent marker-weight perturbations."
    ))
  }
  shrinkage <- as.integer(shrinkage)
  if (length(shrinkage) != 1L || is.na(shrinkage) || !(shrinkage %in% 0:2)) {
    cli::cli_abort(c(
      x = "{.arg shrinkage} must be one of 0, 1, or 2.",
      i = "The mappings are 0 = Student-t(6), 1 = Laplace, and 2 = Normal."
    ))
  }
  compat_args <- list(...)
  unknown_args <- setdiff(names(compat_args), c("n_obs_per_marker_per_id", "n_t", ""))
  if (length(unknown_args) > 0L) {
    cli::cli_abort(c(
      x = "Unknown argument{?s}: {.field {unknown_args}}.",
      i = "Only compatibility arguments {.field n_obs_per_marker_per_id} and {.field n_t} are accepted via {.arg ...}."
    ))
  }
  obs_time_noise_sd <- as.numeric(obs_time_noise_sd %||% 0)
  left_truncation_max <- as.numeric(left_truncation_max %||% 0)
  if (length(left_truncation_max) != 1L || !is.finite(left_truncation_max) || left_truncation_max < 0) {
    cli::cli_abort(c(
      x = "{.arg left_truncation_max} must be one non-negative finite number.",
      i = "Use 0 to disable delayed entry in simulated event data."
    ))
  }
  if (!is.logical(truncate_longitudinal_before_entry) || length(truncate_longitudinal_before_entry) != 1L || is.na(truncate_longitudinal_before_entry)) {
    cli::cli_abort(c(
      x = "{.arg truncate_longitudinal_before_entry} must be TRUE/FALSE.",
      i = "Set TRUE to drop longitudinal rows before delayed-entry times."
    ))
  }
  if (length(obs_time_noise_sd) != 1L || !is.finite(obs_time_noise_sd) || obs_time_noise_sd < 0) {
    cli::cli_abort(c(
      x = "{.arg obs_time_noise_sd} must be one non-negative finite number.",
      i = "Use 0 to disable observation-time jitter."
    ))
  }
  if (!is.logical(censor_longitudinal_after_event) || length(censor_longitudinal_after_event) != 1L || is.na(censor_longitudinal_after_event)) {
    cli::cli_abort(c(
      x = "{.arg censor_longitudinal_after_event} must be TRUE/FALSE.",
      i = "Set TRUE to truncate longitudinal measurements at the event time, or FALSE to keep them through {.arg time_cens}."
    ))
  }
  times_obs <- as.numeric(times_obs %||% numeric(0))
  times_obs <- times_obs[is.finite(times_obs)]
  if (length(times_obs) == 0L) {
    cli::cli_abort(c(
      x = "{.arg times_obs} must contain at least one finite scheduled observation time.",
      i = "The number of longitudinal observations is now inferred directly from {.arg times_obs}."
    ))
  }
  times_obs <- sort(times_obs)

  .sim_surv_lhs_arity <- function(formula_event) {
    lhs <- formula_event[[2]]
    if (!is.call(lhs)) return(NA_integer_)
    head_expr <- lhs[[1]]
    is_surv <- identical(as.character(head_expr), "Surv") ||
      (is.call(head_expr) && identical(as.character(head_expr[[1]]), "::") && identical(as.character(head_expr[[3]]), "Surv"))
    if (!is_surv) return(NA_integer_)
    max(0L, length(lhs) - 1L)
  }
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

  # Reuse the package-neutral interpreter module.  Binding the functions into
  # this simulation closure also makes the complete evaluator available when
  # the closure is serialised to a mirai worker.
  .sim_eval_bytecode_scalar <- eval_bytecode_scalar
  .sim_eval_bytecode_vector <- eval_bytecode_vector

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
    tv_specs <- list()

    .sim_draw_time_varyring <- function(fun, n_step = NULL, steps = NULL, fun_args, n_subjects, time_horizon) {
      steps <- steps %||% NULL
      steps_clean <- NULL
      if (!is.null(steps)) {
        steps_num <- as.numeric(steps)
        if (length(steps_num) == 0L || any(!is.finite(steps_num)) || any(steps_num <= 0) || any(abs(steps_num - round(steps_num)) > 1e-8)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: {.arg steps} must be a vector of positive integers.",
            i = "Example: {.code steps = c(1, 3, 5)}."
          ))
        }
        steps_clean <- sort(unique(as.integer(round(steps_num))))
      } else {
        n_step <- as.integer(n_step)
        if (!(length(n_step) %in% c(1L, 2L)) || any(!is.finite(n_step)) || any(n_step < 1L)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring} requires {.arg n_step} with length 1 or 2 and positive integers when {.arg steps} is NULL.",
            i = "Use {.code n_step = 4} or {.code n_step = c(3, 7)}."
          ))
        }
      }

      if (is.null(steps_clean)) {
        if (length(n_step) == 1L) {
          n_steps_by_id <- rep.int(n_step, n_subjects)
        } else {
          n_min <- min(n_step)
          n_max <- max(n_step)
          n_steps_by_id <- sample.int(n_max - n_min + 1L, size = n_subjects, replace = TRUE) + n_min - 1L
        }
      } else {
        n_steps_by_id <- rep.int(length(steps_clean) + 1L, n_subjects)
      }

      out_specs <- vector("list", n_subjects)
      for (ii in seq_len(n_subjects)) {
        n_segments <- max(1L, as.integer(n_steps_by_id[ii]))
        if (!is.null(steps_clean)) {
          steps_i <- steps_clean[steps_clean < time_horizon]
          breaks <- sort(unique(c(0, steps_i, time_horizon)))
          if (length(breaks) < 2L) {
            breaks <- c(0, time_horizon)
          }
          n_segments <- length(breaks) - 1L
        } else {
          if (n_segments <= 1L) {
            breaks <- c(0, time_horizon)
          } else {
            inner <- sort(stats::runif(n_segments - 1L, min = 0, max = time_horizon))
            breaks <- c(0, inner, time_horizon)
          }
        }
        n_draw <- length(breaks) - 1L
        vals <- do.call(fun, c(list(n = n_draw), fun_args))
        if (length(vals) == 1L) vals <- rep(vals, n_draw)
        if (length(vals) != n_draw) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: generator returned length {length(vals)} for {.val {n_draw}} segments.",
            i = "Ensure {.arg fun} returns one value per segment (or a scalar)."
          ))
        }
        out_specs[[ii]] <- list(
          breaks = as.numeric(breaks),
          values = as.numeric(vals),
          steps = if (is.null(steps_clean)) integer(0) else as.integer(steps_clean)
        )
      }

      out_specs
    }

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

      is_tv_call <- is.call(rhs_expr) && identical(as.character(rhs_expr[[1]]), "time_varyring")
      if (is_tv_call) {
        rhs_args <- as.list(rhs_expr)[-1]
        rhs_arg_names <- names(rhs_args) %||% rep("", length(rhs_args))

        idx_fun <- which(rhs_arg_names == "fun")[1]
        if (is.na(idx_fun)) idx_fun <- 1L
        fn <- eval(rhs_args[[idx_fun]], envir = eval_env)
        if (!is.function(fn)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: {.arg fun} must evaluate to a function.",
            i = "Use generators such as {.fn rnorm}, {.fn rbinom}, or a custom function."
          ))
        }

        idx_steps <- which(rhs_arg_names == "steps")[1]
        idx_n_step <- which(rhs_arg_names == "n_step")[1]
        if (is.na(idx_n_step)) {
          pos_candidates <- setdiff(seq_along(rhs_args), c(idx_fun, idx_steps))
          if (length(pos_candidates) > 0L) idx_n_step <- pos_candidates[1]
        }

        steps <- if (!is.na(idx_steps)) eval(rhs_args[[idx_steps]], envir = eval_env) else NULL
        n_step <- if (!is.na(idx_n_step)) eval(rhs_args[[idx_n_step]], envir = eval_env) else NULL
        if (is.null(steps) && is.null(n_step)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring} requires either {.arg n_step} or {.arg steps}.",
            i = "Example: {.code x ~ time_varyring(rnorm, c(3, 6), mean = 0, sd = 1)} or {.code x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1)}."
          ))
        }

        consumed <- unique(c(idx_fun, idx_n_step, idx_steps))
        consumed <- consumed[is.finite(consumed) & consumed >= 1L]
        extra_idx <- setdiff(seq_along(rhs_args), consumed)
        extra_args <- if (length(extra_idx) > 0L) lapply(rhs_args[extra_idx], eval, envir = eval_env) else list()

        tv_specs[[lhs_name]] <- .sim_draw_time_varyring(
          fun = fn,
          n_step = n_step,
          steps = steps,
          fun_args = extra_args,
          n_subjects = n_subjects,
          time_horizon = time_cens
        )
        sampled <- vapply(tv_specs[[lhs_name]], function(spec) spec$values[1], numeric(1))
      } else {
        sampled <- eval(rhs_expr, envir = eval_env)
      }

      if (length(sampled) == 1) sampled <- rep(sampled, n_subjects)
      if (length(sampled) != n_subjects) {
        cli::cli_abort(c(
          x = "Covariate formula for '{lhs_name}' returned length {length(sampled)}.",
          i = "Expected length {.val {n_subjects}} or scalar."
        ))
      }
      out[[lhs_name]] <- sampled
    }
    attr(out, "time_varyring_specs") <- tv_specs
    out
  }

  .sim_draw_re_block <- function(n_group, K, cfg, label, resolved = NULL) {
    if (K <= 0) return(matrix(0.0, n_group, 0))
    resolved <- resolved %||% .sim_resolve_re_block_cfg(cfg, K, label)
    sd_vec <- resolved$sd
    corr <- resolved$corr
    Sigma <- diag(as.numeric(sd_vec), K, K) %*% corr %*% diag(as.numeric(sd_vec), K, K)
    .sim_draw_mvn(n_group, Sigma)
  }

  .sim_draw_mvn <- function(n_group, Sigma) {
    k_dim <- ncol(Sigma)
    if (k_dim <= 0L || n_group <= 0L) {
      return(matrix(0.0, n_group, k_dim))
    }
    R <- chol(Sigma)
    Z <- matrix(stats::rnorm(n_group * k_dim), n_group, k_dim)
    Z %*% R
  }

  .sim_rhs_matrix <- function(rhs_list, data) {
    if (length(rhs_list) == 0) {
      return(matrix(0.0, nrow(data), 0))
    }
    mats <- lapply(rhs_list, function(rhs) .mm(rhs, data))
    do.call(cbind, mats)
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
      bytecode = inv_link_bc$bytecode,
      const_data = inv_link_bc$const_data %||% numeric(0)
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

  .sim_validate_fixed_transform_shifts <- function(transforms) {
    if (is.null(transforms) || !length(transforms)) {
      return(transforms)
    }
    has_fit_only_iota <- vapply(transforms, function(spec) {
      isTRUE((spec$n_iota_intercept %||% 0L) > 0L) ||
        isTRUE((spec$n_iota_slope %||% 0L) > 0L) ||
        isTRUE(spec$estimate_iota_intercept) ||
        isTRUE(spec$estimate_iota_slope)
    }, logical(1))
    if (any(has_fit_only_iota)) {
      bad_terms <- names(transforms)[has_fit_only_iota]
      cli::cli_abort(c(
        x = "Fit-only functional transform affine shifts are not available in {.fn simulate_joinme}.",
        i = "Remove {.arg intercept}/{.arg slope} from {.fn joinme_tf} for simulation, or encode fixed shifts directly in the transform expression for: {.val {paste(bad_terms, collapse = ', ')}}."
      ))
    }
    transforms
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
          bytecode = bc$bytecode,
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
            i = "Stan-estimated penalised splines without y are supported in {.fn JoiNMe}, not in simulation."
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
    tag_truth <- function(fun, type, parameters = list(), formula = NULL, coefficients = NULL) {
      attr(fun, "basehaz_truth") <- list(
        type = type,
        parameters = parameters,
        formula = formula,
        coefficients = coefficients
      )
      fun
    }
    if (!is.null(formula_basehaz)) {
      f_bh <- stats::as.formula(formula_basehaz)
      t_ref <- unique(c(0, time_cens / 2, time_cens))
      X_ref <- stats::model.matrix(f_bh, data = data.frame(time = t_ref))
      coef_bh <- .sim_align_coef(
        colnames(X_ref),
        user_coef = beta_basehaz,
        sd_default = 0.2,
        intercept_default = -2.0
      )

      fun <- function(t) {
        df_t <- data.frame(time = as.numeric(t))
        X_bh <- stats::model.matrix(f_bh, data = df_t)
        lp <- as.numeric(X_bh %*% coef_bh)
        lp <- pmin(lp, log(.Machine$double.xmax) - 2)
        as.numeric(exp(lp))
      }
      return(tag_truth(fun, "formula", formula = f_bh, coefficients = coef_bh))
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
      rate <- as.numeric(spec$rate %||% spec$lambda %||% runif(1, min = 0.01, max = 0.5))
      if (isTRUE(spec$log)) {
        rate <- exp(rate)
      }
      fun <- function(t) rep(rate, length(t))
      return(tag_truth(
        fun,
        "constant",
        parameters = list(rate = rate),
        formula = ~1,
        coefficients = stats::setNames(log(rate), "(Intercept)")
      ))
    }

    if (mode == "linear") {
      intercept <- as.numeric(spec$intercept %||% -2.2)
      slope <- as.numeric(spec$slope %||% 0.25)
      fun <- function(t) {
        t <- as.numeric(t)
        # softplus(intercept + slope * pmax(t, 0))
        exp(intercept + slope * pmax(t, 0))
      }
      return(tag_truth(
        fun,
        "linear",
        parameters = list(intercept = intercept, slope = slope),
        formula = ~1 + time,
        coefficients = c("(Intercept)" = intercept, time = slope)
      ))
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
      fun <- function(t) {
        t <- pmax(as.numeric(t), 0)
        idx <- findInterval(t, vec = breaks, rightmost.closed = TRUE) + 1L
        rates[idx]
      }
      return(tag_truth(fun, "piecewise", parameters = list(breaks = breaks, rates = rates)))
    }

    if (mode == "weibull") {
      shape <- as.numeric(spec$shape %||% 1.4)
      scale <- as.numeric(spec$scale %||% 6.0)
      return(tag_truth(
        weibull_h0(shape = shape, scale = scale),
        "weibull",
        parameters = list(shape = shape, scale = scale)
      ))
    }

    if (mode %in% c("spline", "bs", "ns")) {
      basis_type <- tolower(as.character(spec$basis %||% if (mode == "ns") "ns" else "bs")[1])
      knots <- as.numeric(spec$knots %||% stats::quantile(seq(0, time_cens, length.out = 100), probs = c(0.25, 0.5, 0.75)))
      degree <- as.integer(spec$degree %||% 3L)
      coef <- spec$coef %||% beta_basehaz
      intercept <- as.logical(spec$intercept %||% TRUE)

      build_basis <- function(t) {
        if (basis_type == "ns") {
          splines::ns(t, knots = knots, Boundary.knots = c(0, time_cens), intercept = intercept)
        } else {
          splines::bs(t, knots = knots, Boundary.knots = c(0, time_cens), degree = degree, intercept = intercept)
        }
      }

      # Resolve spline coefficients once so baseline hazard is deterministic over
      # repeated evaluations during integration/root-finding.
      t_ref <- unique(c(0, time_cens / 2, time_cens))
      B_ref <- as.matrix(build_basis(t_ref))
      coef_vec <- .sim_align_coef(
        colnames(B_ref),
        user_coef = coef,
        sd_default = 0.25,
        intercept_default = -2.2
      )

      fun <- function(t) {
        t <- pmax(as.numeric(t), 0)
        B <- as.matrix(build_basis(t))
        lp <- as.numeric(B %*% coef_vec)
        lp <- pmin(lp, log(.Machine$double.xmax) - 2)
        as.numeric(exp(lp))
      }
      return(tag_truth(
        fun,
        basis_type,
        parameters = list(knots = knots, degree = degree, intercept = intercept),
        coefficients = coef_vec
      ))
    }

    cli::cli_abort(c(
      x = "Unsupported baseline hazard mode: {.val {mode}}.",
      i = "Use one of: constant, linear, piecewise, weibull, spline."
    ))
  }

  # ---- Step 1: Define marker/family dimensions and parse model structure
  if (.is_single_family_spec(families)) {
    n_markers_requested <- length(marker_levels %||% "m1")
    families <- rep(list(families), n_markers_requested)
  } else if (
    length(families) == 1L &&
      !is.null(marker_levels) &&
      length(marker_levels) > 1L
  ) {
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
  # ---- Step 1a: Resolve marker-family mappings and inverse-link bytecode
  family_spec <- .sim_resolve_family_specs(families, D)
  family_codes <- family_spec$family_codes
  link_names <- family_spec$link_names
  inv_link_specs <- family_spec$inv_link_specs
  use_tau_fixed <- family_spec$use_tau_fixed
  tau_fixed <- family_spec$tau_fixed
  family_names <- vapply(family_codes, .family_code_to_name, character(1))
  skew_laplace_marker <- family_codes == 9L

  # Resolve distributional formulas before generating any random quantities.
  # This permits immediate rejection of a tau regression that no simulated
  # marker can use, avoiding unnecessary longitudinal and event-time sampling.
  dist_formulas <- .normalize_formula_dist(formulaDist)
  family_names_present <- vapply(
    sort(unique(family_codes)),
    .family_code_to_name,
    character(1)
  )
  .validate_dist_formula_scopes(dist_formulas, family_names_present)
  if (
    !is.null(dist_formulas$tau) &&
      any(skew_laplace_marker) &&
      all(use_tau_fixed[skew_laplace_marker] == 1L)
  ) {
    cli::cli_abort(c(
      x = "The {.code tau ~ ...} distributional regression has no simulated skew-Laplace marker.",
      i = "Remove the tau regression or omit {.arg tau} from at least one skew-Laplace family specification."
    ))
  }

  # ---- Step 1b: Parse longitudinal random-effect structure from formulaLong
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
  indep_flags <- .get_re_independence(formulaLong, marker_var = marker_var, id_var = id_var)
  id_idx <- which(grp_names == id_var)
  id_rhs_list <- if (length(id_idx) > 0) .bar_terms_to_rhs_list(bars[id_idx]) else list()

  nested_terms <- .extract_nested_marker_terms(formulaLong, marker_var = marker_var, id_var = id_var)
  mk_rhs_list <- nested_terms$mk_rhs_list
  idm_rhs_list <- nested_terms$idm_rhs_list

  # ---- Step 2: Generate exogenous event-level covariates and event design
  dataEvent <- data.frame(id = seq_len(n_id), stringsAsFactors = FALSE)
  dataEvent <- .sim_eval_covariate_formulas(n_subjects = n_id, formulas = covariate_formulas, base_df = dataEvent)
  tv_cov_specs <- attr(dataEvent, "time_varyring_specs") %||% list()
  attr(dataEvent, "time_varyring_specs") <- NULL

  surv_lhs_arity <- .sim_surv_lhs_arity(formulaEvent)
  split_required <- (isTRUE(surv_lhs_arity >= 3L) || length(tv_cov_specs) > 0L || left_truncation_max > 0)
  event_layout <- if (split_required) "split" else "subject"

  if (identical(event_layout, "split") && !isTRUE(surv_lhs_arity >= 3L)) {
    cli::cli_abort(c(
      x = "Automatic split event rows require counting-process survival syntax in {.arg formulaEvent}.",
      i = "Use {.code survival::Surv(time_start, time_stop, event) ~ ...}."
    ))
  }
  names(dataEvent)[names(dataEvent) == "id"] <- id_var

  # Covariance-regression design for the subject-specific marker-by-id covariance.
  formulaVCov <- .get_vcov_formula(
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

  # ---- Step 2a: Build prototype matrices for deterministic coefficient alignment
  # Time-dependent model matrices are defined on the scaled [0, 1] domain so
  # spline bases follow the same construction path as joinme_standata().
  time_scale_internal <- max(c(time_cens, times_obs), na.rm = TRUE)
  if (!is.finite(time_scale_internal) || time_scale_internal <= 0) {
    time_scale_internal <- 1.0
  }

  prototype <- dataEvent[rep(1, D), , drop = FALSE]
  prototype[[marker_var]] <- factor(marker_levels, levels = marker_levels)
  prototype[[time_var]] <- rep(0.5, D)

  fixed_template <- .make_model_matrix_template(
    fixed_rhs,
    prototype,
    boundary_var = time_var,
    boundary_values = c(0, 1)
  )
  id_templates <- lapply(id_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype, boundary_var = time_var, boundary_values = c(0, 1))
  })
  mk_templates <- lapply(mk_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype, boundary_var = time_var, boundary_values = c(0, 1))
  })
  idm_templates <- lapply(idm_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype, boundary_var = time_var, boundary_values = c(0, 1))
  })

  X_proto <- .mm(fixed_template, prototype)
  Z_id_proto <- .sim_rhs_matrix(id_templates, prototype)
  Z_mk_proto <- .sim_rhs_matrix(mk_templates, prototype)
  Z_idm_proto <- .sim_rhs_matrix(idm_templates, prototype)

  time_meta <- .make_time_index_metadata(
    time_var = time_var,
    fixed_design = fixed_template,
    id_design = id_templates,
    marker_design = mk_templates,
    idm_design = idm_templates
  )
  idx_time_beta <- time_meta$idx_time_beta
  idx_time_uid <- time_meta$idx_time_uid
  idx_time_vmk <- time_meta$idx_time_vmk
  idx_time_idm <- time_meta$idx_time_idm

  beta_long <- .sim_align_coef(colnames(X_proto), beta_long, sd_default = 0.35, intercept_default = 1.0)
  beta_long_internal <- beta_long
  if (length(idx_time_beta) > 0L) {
    beta_long_internal[idx_time_beta] <- beta_long_internal[idx_time_beta] * time_scale_internal
  }

  assoc_from_formula <- .sim_assoc_from_formula(formulaAssoc)
  assoc_effective <- if (!is.null(assoc_from_formula)) assoc_from_formula else assoc
  assoc_effective <- unique(assoc_effective)
  if (length(assoc_effective) == 0) assoc_effective <- "cv_total"

  marker_weight_spec <- .get_marker_weight_structure(
    marker_weights = marker_weights,
    marker_levels = marker_levels,
    active_terms = .active_weighted_assoc_terms(assoc_effective),
    shared_marker_weights = shared_marker_weights,
    estimate_marker_weights = !isTRUE(fixed_marker_weights),
    context = "simulate_joinme()"
  )
  marker_weights_base_by_term <- marker_weight_spec$base_by_term
  z_marker_weight_sets <- matrix(
    0,
    nrow = marker_weight_spec$n_sets,
    ncol = D,
    dimnames = dimnames(marker_weight_spec$base_matrix)
  )
  if (!isTRUE(fixed_marker_weights) && length(marker_weight_spec$active_terms) > 0L) {
    z_marker_weight_sets[] <- .sim_draw_standard_shrinkage(
      marker_weight_spec$n_sets * D,
      shrinkage = shrinkage
    )
  }
  marker_weights_eff <- marker_weights_base_by_term
  marker_weights_latent_by_term <- marker_weights_base_by_term
  for (term_key in names(marker_weights_latent_by_term)) {
    marker_weights_latent_by_term[[term_key]][] <- 0
    set_pos <- marker_weight_spec$set_index[[term_key]]
    if (set_pos > 0L) {
      marker_weights_latent_by_term[[term_key]] <- as.numeric(z_marker_weight_sets[set_pos, ])
      marker_weights_eff[[term_key]] <-
        as.numeric(marker_weights_base_by_term[[term_key]]) +
        marker_weights_latent_by_term[[term_key]]
    }
    names(marker_weights_latent_by_term[[term_key]]) <- marker_levels
    names(marker_weights_eff[[term_key]]) <- marker_levels
  }
  .sim_public_marker_weights <- function(by_term) {
    if (isTRUE(marker_weight_spec$shared_marker_weights)) {
      term_key <- if (length(marker_weight_spec$active_terms) > 0L) {
        marker_weight_spec$active_terms[[1L]]
      } else {
        .weighted_assoc_term_keys()[[1L]]
      }
      return(by_term[[term_key]])
    }
    by_term[marker_weight_spec$active_terms]
  }
  marker_weights_base <- .sim_public_marker_weights(marker_weights_base_by_term)
  marker_weights_latent <- .sim_public_marker_weights(marker_weights_latent_by_term)
  marker_weights_effective <- .sim_public_marker_weights(marker_weights_eff)

  .sim_time_varyring_value <- function(spec_i, t_val) {
    breaks <- as.numeric(spec_i$breaks)
    values <- as.numeric(spec_i$values)
    t_use <- min(max(as.numeric(t_val), breaks[1]), breaks[length(breaks)])
    seg <- findInterval(t_use, breaks, rightmost.closed = TRUE, all.inside = TRUE)
    seg <- min(max(seg, 1L), length(values))
    values[seg]
  }

  .sim_event_row_at <- function(i, t_val) {
    row <- dataEvent[i, , drop = FALSE]
    if (length(tv_cov_specs) > 0L) {
      for (nm in names(tv_cov_specs)) {
        row[[nm]] <- .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], t_val)
      }
    }
    row
  }

  event_template <- do.call(rbind, lapply(seq_len(n_id), function(i) .sim_event_row_at(i, 0)))
  W_event <- .mm_event(formulaEvent, event_template)

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

  # ---- Step 3: Draw hierarchical random effects and covariance-regression factors
  K_id <- ncol(Z_id_proto)
  K_mk <- ncol(Z_mk_proto)
  K_idm <- ncol(Z_idm_proto)

  .validate_assoc_cov_structure(
    assoc = assoc_effective,
    q_idm = K_idm,
    diagonal_only = indep_flags$indep_idmarker_cov,
    context = "simulate_joinme()"
  )

  re_id_effective <- .sim_resolve_re_block_cfg(re_params$id %||% list(), K_id, "id")
  re_marker_effective <- .sim_resolve_re_block_cfg(re_params$marker %||% list(), K_mk, "marker")
  if (K_id > 1L && as.integer(indep_flags$indep_id_re %||% 0L) == 1L) {
    id_corr_offdiag <- re_id_effective$corr
    diag(id_corr_offdiag) <- 0
    if (any(abs(id_corr_offdiag) > 1e-8)) {
      cli::cli_warn(c(
        x = "Top-level {.code || id} requests independent subject-level random effects in {.fn simulate_joinme}.",
        i = "Ignoring nonzero off-diagonal entries in {.arg re_params$id$corr}."
      ))
    }
    re_id_effective$corr <- diag(K_id)
    re_id_effective$cov <- diag(as.numeric(re_id_effective$sd), K_id, K_id)^2
    re_id_effective$Lcorr <- diag(K_id)
  }
  re_id_effective_internal <- re_id_effective
  if (length(idx_time_uid) > 0L) {
    re_id_effective_internal$sd[idx_time_uid] <- re_id_effective_internal$sd[idx_time_uid] * time_scale_internal
    re_id_effective_internal$cov <- diag(as.numeric(re_id_effective_internal$sd), K_id, K_id) %*%
      re_id_effective_internal$corr %*%
      diag(as.numeric(re_id_effective_internal$sd), K_id, K_id)
  }
  re_marker_effective_internal <- re_marker_effective
  if (length(idx_time_vmk) > 0L) {
    re_marker_effective_internal$sd[idx_time_vmk] <- re_marker_effective_internal$sd[idx_time_vmk] * time_scale_internal
    re_marker_effective_internal$cov <- diag(as.numeric(re_marker_effective_internal$sd), K_mk, K_mk) %*%
      re_marker_effective_internal$corr %*%
      diag(as.numeric(re_marker_effective_internal$sd), K_mk, K_mk)
  }
  # ---- Step 3a: prepare the optional latent-progress distribution
  # The mixture builder is invoked only after the random-effect formulae have
  # determined every available dimension. This is the earliest point at which
  # `class_dimensions` can be checked faithfully. The builder itself reuses
  # the fitting coordinate-layout machinery, so simulation and estimation
  # cannot silently disagree about which standardised coordinate an index
  # denotes.
  mixture_configuration <- .sim_prepare_mixture(
    specification = .mixture_specification,
    data_event = dataEvent,
    marker_prototype = prototype,
    dimensions = c(
      subject = K_id,
      marker = K_mk,
      covariance_basis = K_idm
    ),
    labels = list(
      subject = colnames(Z_id_proto) %||% character(0),
      marker = colnames(Z_mk_proto) %||% character(0),
      covariance = colnames(Z_idm_proto) %||% character(0)
    ),
    covariance_is_diagonal =
      as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L,
    shrinkage = shrinkage,
    id_variable = id_var,
    marker_variable = marker_var
  ) # checked component layout, probabilities, allocations, locations and scales

  # Draw the standardised individual effects first, replace only the selected
  # coordinates conditional on the subject's common class, and finally apply
  # the ordinary random-effect covariance factor. This is precisely the
  # `u_i = L_u z_{u,i}` ordering used in Stan.
  z_id <- if (K_id > 0L) {
    matrix(stats::rnorm(n_id * K_id), nrow = n_id, ncol = K_id)
  } else {
    matrix(0, nrow = n_id, ncol = 0L)
  } # standardised individual random effects before covariance scaling
  z_id <- .sim_apply_mixture_to_latent(
    latent_matrix = z_id,
    level = "subject",
    allocation = mixture_configuration$allocation$subject %||% integer(0),
    mixture = mixture_configuration,
    shrinkage = shrinkage
  ) # selected component-conditional individual coordinates
  re_id <- if (K_id > 0L) {
    z_id %*% chol(re_id_effective_internal$cov)
  } else {
    matrix(0, nrow = n_id, ncol = 0L)
  } # realised individual effects on the internally scaled time basis

  # Repeat the same construction for marker-level standardised effects. A
  # combined marker and marker-weight mixture uses the already drawn common
  # marker allocation rather than sampling a second label.
  z_marker <- if (K_mk > 0L) {
    matrix(stats::rnorm(D * K_mk), nrow = D, ncol = K_mk)
  } else {
    matrix(0, nrow = D, ncol = 0L)
  } # standardised marker random effects before covariance scaling
  z_marker <- .sim_apply_mixture_to_latent(
    latent_matrix = z_marker,
    level = "marker",
    allocation = mixture_configuration$allocation$marker %||% integer(0),
    mixture = mixture_configuration,
    shrinkage = shrinkage
  ) # selected component-conditional marker coordinates
  re_marker <- if (K_mk > 0L) {
    z_marker %*% chol(re_marker_effective_internal$cov)
  } else {
    matrix(0, nrow = D, ncol = 0L)
  } # realised marker effects on the internally scaled time basis

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
    .sim_draw_mvn(n_id * D, diag(K_idm))
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
    .cov_lp_to_chol(
      lp_vec = lp_vec,
      q_idm = q_idm,
      idx_row = idx_row,
      idx_col = idx_col,
      diag_link = diag_link
    )
  }

  .sim_cov_matrix_to_alpha <- function(mat, idx_row, idx_col, diag_link) {
    .cov_chol_to_lp(
      L_i = mat,
      idx_row = idx_row,
      idx_col = idx_col,
      diag_link = diag_link
    )
  }

  .sim_align_len <- function(x, n, default, arg_name = "parameter") {
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
    if (length(x) < n) {
      cli::cli_abort(c(
        x = "Length mismatch in covariance-regression parameter specification.",
        i = "Expected length {n}, got {length(x)}."
      ))
    }
    if (length(x) > n) {
      cli::cli_warn(c(
        x = "Ignoring extra covariance-regression values in {.arg {arg_name}}.",
        i = "The current marker-by-id covariance structure uses {n} component{?s}, but {length(x)} value{?s} were supplied."
      ))
      return(x[seq_len(n)])
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
  alpha_cov <- .sim_align_len(re_cov_cfg$alpha, M_cov, default = default_alpha_cov, arg_name = "re_params$id_marker_cov$alpha")
  if (!is.null(re_cov_cfg$lambda) && !is.null(dim(re_cov_cfg$lambda))) {
    cli::cli_abort(c(
      x = "{.arg re_params$id_marker_cov$lambda} must be a scalar or numeric vector.",
      i = "A matrix would imply cross-coordinate loadings, which this element-wise covariance regression does not use."
    ))
  }
  lambda_cov <- .sim_align_len(re_cov_cfg$lambda, M_cov, default = default_lambda_cov, arg_name = "re_params$id_marker_cov$lambda")
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
  marker_id_row_scale_internal <- rep(1.0, K_idm)
  if (length(idx_time_idm) > 0L) {
    marker_id_row_scale_internal[idx_time_idm] <- time_scale_internal
  }
  z_cov <- matrix(0.0, nrow = n_id, ncol = max(1L, M_cov))
  lambda_cov_sign <- rep(1, M_cov)
  if (M_cov > 0 && K_idm > 0) {
    z_cov_raw <- matrix(stats::rnorm(n_id * M_cov), nrow = n_id, ncol = M_cov)
    cov_latent <- .sim_canonicalize_cov_latent(lambda_cov, z_cov_raw)
    lambda_cov <- cov_latent$lambda
    z_cov <- cov_latent$z
    lambda_cov_sign <- cov_latent$sign

    # The covariance-regression mixture acts on Stan's canonical `z_L`
    # coordinates, after any sign absorbed from a legacy negative loading.
    # Individual and covariance-regression classes share the same subject
    # allocation held in `mixture_configuration$allocation$subject`.
    covariance_class_type <- intersect(
      c("corr", "vcov"),
      mixture_configuration$class_type %||% character(0)
    ) # public covariance representation selected for the shared latent block
    z_cov <- .sim_apply_mixture_to_latent(
      latent_matrix = z_cov,
      level = if (length(covariance_class_type) == 1L) {
        covariance_class_type[[1L]]
      } else {
        ""
      },
      allocation =
        mixture_configuration$allocation$subject %||% integer(0),
      mixture = mixture_configuration,
      shrinkage = shrinkage
    ) # selected component-conditional covariance-regression latents

    lp_cov <- matrix(alpha_cov, nrow = n_id, ncol = M_cov, byrow = TRUE)
    if (K_cov > 0) {
      lp_cov <- lp_cov + Xcov %*% t(beta_cov)
    }
    lp_cov <- lp_cov + sweep(z_cov, 2L, lambda_cov, `*`)

    for (i in seq_len(n_id)) {
      L_i[i, , ] <- .sim_cov_lp_to_matrix(
        lp_vec = lp_cov[i, ],
        q_idm = K_idm,
        idx_row = idx_row_cov,
        idx_col = idx_col_cov,
        diag_link = diag_link_cov
      )
    }
  }

  L_i_eff <- L_i
  if (K_idm > 0L && any(marker_id_row_scale_internal != 1)) {
    for (i in seq_len(n_id)) {
      L_i_eff[i, , ] <- sweep(L_i_eff[i, , ], 1L, marker_id_row_scale_internal, `*`)
    }
  }

  re_idm_scaled <- array(0.0, dim = c(n_id, D, K_idm))
  if (K_idm > 0) {
    for (i in seq_len(n_id)) {
      Li <- matrix(L_i_eff[i, , ], K_idm, K_idm)
      re_idm_scaled[i, , ] <- re_idm[i, , ] %*% t(Li)
    }
  }

  rescale_re_matrix_to_public <- function(mat, idx, scale_factor) {
    if (!is.matrix(mat) || ncol(mat) == 0L) return(mat)
    idx <- as.integer(idx %||% integer(0))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= ncol(mat)]
    if (!length(idx) || !is.finite(scale_factor) || scale_factor <= 0 || abs(scale_factor - 1) < 1e-12) {
      return(mat)
    }
    mat[, idx] <- mat[, idx, drop = FALSE] / scale_factor
    mat
  }

  rescale_re_array_to_public <- function(arr, idx, scale_factor) {
    if (is.null(arr) || length(dim(arr)) != 3L || dim(arr)[3] == 0L) return(arr)
    idx <- as.integer(idx %||% integer(0))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= dim(arr)[3]]
    if (!length(idx) || !is.finite(scale_factor) || scale_factor <= 0 || abs(scale_factor - 1) < 1e-12) {
      return(arr)
    }
    arr[, , idx] <- arr[, , idx, drop = FALSE] / scale_factor
    arr
  }

  label_re_matrix <- function(mat, row_labels, col_labels) {
    if (!is.matrix(mat)) {
      return(mat)
    }
    row_labels <- as.character(row_labels %||% character(0))
    col_labels <- as.character(col_labels %||% character(0))
    if (nrow(mat) == length(row_labels)) {
      rownames(mat) <- row_labels
    }
    if (ncol(mat) == length(col_labels)) {
      colnames(mat) <- col_labels
    }
    mat
  }

  label_re_array <- function(arr, id_labels, marker_labels, term_labels) {
    if (is.null(arr) || length(dim(arr)) != 3L) {
      return(arr)
    }

    dim_names <- dimnames(arr)
    if (is.null(dim_names)) {
      dim_names <- vector("list", 3L)
    }

    id_labels <- as.character(id_labels %||% character(0))
    marker_labels <- as.character(marker_labels %||% character(0))
    term_labels <- as.character(term_labels %||% character(0))

    if (dim(arr)[1] == length(id_labels)) {
      dim_names[[1]] <- id_labels
    }
    if (dim(arr)[2] == length(marker_labels)) {
      dim_names[[2]] <- marker_labels
    }
    if (dim(arr)[3] == length(term_labels)) {
      dim_names[[3]] <- term_labels
    }

    dimnames(arr) <- dim_names
    arr
  }

  re_id_public <- rescale_re_matrix_to_public(re_id, idx_time_uid, time_scale_internal)
  re_marker_public <- rescale_re_matrix_to_public(re_marker, idx_time_vmk, time_scale_internal)
  re_idm_public <- rescale_re_array_to_public(re_idm_scaled, idx_time_idm, time_scale_internal)

  id_labels_public <- as.character(dataEvent[[id_var]] %||% seq_len(n_id))
  id_term_labels <- colnames(Z_id_proto) %||% paste0("id_re_", seq_len(ncol(re_id_public)))
  marker_term_labels <- colnames(Z_mk_proto) %||% paste0("marker_re_", seq_len(ncol(re_marker_public)))
  idm_term_labels <- colnames(Z_idm_proto) %||% paste0("id_marker_re_", seq_len(dim(re_idm_public)[3]))

  re_id <- label_re_matrix(re_id, id_labels_public, id_term_labels)
  re_id_public <- label_re_matrix(re_id_public, id_labels_public, id_term_labels)
  re_marker <- label_re_matrix(re_marker, marker_levels, marker_term_labels)
  re_marker_public <- label_re_matrix(re_marker_public, marker_levels, marker_term_labels)
  re_idm <- label_re_array(re_idm, id_labels_public, marker_levels, idm_term_labels)
  re_idm_scaled <- label_re_array(re_idm_scaled, id_labels_public, marker_levels, idm_term_labels)
  re_idm_public <- label_re_array(re_idm_public, id_labels_public, marker_levels, idm_term_labels)

  # ---- Step 3a: Build association channels entering the survival linear predictor
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

  transforms <- .sim_validate_fixed_transform_shifts(.sim_normalize_transform_list(transforms))
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

  # ---- Step 3b: Define latent longitudinal trajectory evaluators over continuous time
 
  assoc_row_template <- vector("list", n_id)
  marker_factor <- factor(marker_levels, levels = marker_levels)
  for (i in seq_len(n_id)) {
    row_df_i <- dataEvent[rep(i, D), , drop = FALSE]
    row_df_i[[marker_var]] <- marker_factor
    assoc_row_template[[i]] <- row_df_i
  }

  eta_components_all_markers <- function(i, t) {
    row_df <- assoc_row_template[[i]]
    row_df[[time_var]] <- t / time_scale_internal

    x_fix <- .mm(fixed_template, row_df)
    z_id <- .sim_rhs_matrix(id_templates, row_df)
    z_mk <- .sim_rhs_matrix(mk_templates, row_df)
    z_idm <- .sim_rhs_matrix(idm_templates, row_df)

    fixed_part <- if (ncol(x_fix) > 0) as.numeric(x_fix %*% beta_long_internal) else rep(0, D)
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
    # Build subject-level raw correlation-factor features in the same lower-
    # triangular ordering used by Stan: (2,1), (3,1), (3,2), ... .
    #
    # Important: this returns the off-diagonal entries of the subject-specific
    # Cholesky-correlation factor K_i, not pairwise correlations from Sigma_i.
    if (K_idm < 2) return(numeric(0))
    Li <- matrix(L_i[i, , ], nrow = K_idm, ncol = K_idm)
    if (!all(is.finite(Li))) return(rep(0.0, M_corr))
    .assoc_corr_features_from_chol(Li)
  }

  .sim_vcov_features <- function(i) {
    if (K_idm < 1) return(numeric(0))
    Li <- matrix(L_i_eff[i, , ], nrow = K_idm, ncol = K_idm)
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

    w_mean <- function(term_key, x) {
      weights_term <- as.numeric(marker_weights_eff[[term_key]] %||% rep(1, D))
      sum(weights_term * x) / D
    }

    # Step B: aggregate raw CV summaries (mean, marker, total).
    cv_mean_raw <- mean(now$mu_mean)
    cv_marker_raw <- w_mean("cv_marker", now$mu_marker)
    cv_total_raw <- w_mean("cv_total", now$mu_total)

    cv_mean <- tf_funs$cv_mean(cv_mean_raw)
    cv_marker <- if (has_tf_cv_marker) w_mean("cv_marker", tf_funs$cv_marker(now$mu_marker)) else cv_marker_raw
    cv_total <- if (has_tf_cv_total) w_mean("cv_total", tf_funs$cv_total(now$mu_total)) else cv_total_raw

    # Step C: compute mean slope via finite difference.
    cs_mean_raw <- (mean(eps$mu_mean) - cv_mean_raw) / eps_cs

    # Step D: marker-level slopes for marker and total components.
    # We keep these at marker resolution so CS transforms are applied per marker
    # before weighted averaging, mirroring Stan semantics.
    cs_marker_raw_by_marker <- (eps$mu_marker - now$mu_marker) / eps_cs
    cs_total_raw_by_marker <- (eps$mu_total - now$mu_total) / eps_cs

    # Step E: CS aggregation semantics (aligned with CV aggregation semantics):
    # transform first at marker level, then weighted-average across markers.
    cs_marker <- w_mean("cs_marker", tf_funs$cs_marker(cs_marker_raw_by_marker))
    cs_total <- w_mean("cs_total", tf_funs$cs_total(cs_total_raw_by_marker))

    # Step F: keep weighted raw summaries for debugging/inspection helpers.
    cs_marker_raw <- w_mean("cs_marker", cs_marker_raw_by_marker)
    cs_total_raw <- w_mean("cs_total", cs_total_raw_by_marker)

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
  baseline_hazard_truth <- attr(h0_fn, "basehaz_truth", exact = TRUE)
  if (is.null(baseline_hazard_truth)) {
    baseline_hazard_truth <- list(
      type = "custom_function",
      parameters = list(),
      formula = NULL,
      coefficients = NULL
    )
  }

  .eta_event_it <- function(i, t) {
    row <- .sim_event_row_at(i, t)
    W_row <- .mm_event(formulaEvent, row)
    w <- rep(0, length(beta_event))
    names(w) <- names(beta_event)
    common <- intersect(colnames(W_row), names(beta_event))
    if (length(common) > 0L) {
      w[common] <- as.numeric(W_row[1, common, drop = TRUE])
    }
    sum(w * beta_event)
  }

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
      return(vapply(t, function(tt) h0_fn(tt) * exp(.eta_event_it(i, tt) + .assoc_lp_cached(i, tt)), numeric(1)))
    }
    h0_fn(t) * exp(.eta_event_it(i, t) + .assoc_lp_cached(i, t))
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

    gk_integral <- function() {
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

    value <- if (identical(method, "integrate")) {
      integrate_ctrl <- integration_control
      integrate_ctrl$method <- NULL
      integrate_ctrl$max_refine <- NULL
      out <- tryCatch(
        do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integrate_ctrl)),
        error = function(e) NULL
      )
      if (is.null(out)) gk_integral() else as.numeric(out$value)
    } else {
      gk_integral()
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

  # ---- Step 4: Draw event/censoring times from the implied hazard model
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

  # ---- Step 5a: Sample delayed-entry times once per subject.
  # Entry is constrained to [0, stop_i) so each subject remains at risk after
  # entry and before the observed event/censoring stop time.
  entry_time_by_id <- vapply(seq_len(n_id), function(i) {
    stop_i <- as.numeric(dataEvent[[event_time_var]][i])
    entry_upper <- min(left_truncation_max, max(0, stop_i - 1e-8))
    if (entry_upper > 0) stats::runif(1, min = 0, max = entry_upper) else 0
  }, numeric(1))

  # ---- Step 5b: Build longitudinal observation schedule conditional on
  # event/censoring and optional delayed-entry truncation.
  cov_names <- setdiff(colnames(dataEvent), c(id_var, event_time_var, event_var))

  .sim_obs_rows_for_id <- function(i) {
    entry_i <- as.numeric(entry_time_by_id[i])
    obs_upper <- max(1e-8, min(dataEvent[[event_time_var]][i], time_cens))
    rows_i <- vector("list", D)
    for (d in seq_len(D)) {
      t_obs <- times_obs

      if (obs_time_noise_sd > 0) {
        noise <- stats::rnorm(length(t_obs), mean = 0, sd = obs_time_noise_sd)
        t_obs <- sort(pmin(time_cens, pmax(0, t_obs + noise)))
      }
      if (isTRUE(censor_longitudinal_after_event)) {
        t_obs <- t_obs[t_obs <= obs_upper]
      }
      if (isTRUE(truncate_longitudinal_before_entry)) {
        t_obs <- t_obs[t_obs >= entry_i]
      }
      if (!length(t_obs)) {
        rows_i[[d]] <- NULL
        next
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
    rows_i <- Filter(Negate(is.null), rows_i)
    if (!length(rows_i)) {
      row_df <- data.frame(
        id = dataEvent[[id_var]][i],
        marker = factor(marker_levels[1], levels = marker_levels),
        time = obs_upper,
        stringsAsFactors = FALSE
      )
      names(row_df)[names(row_df) == "id"] <- id_var
      names(row_df)[names(row_df) == "time"] <- time_var
      names(row_df)[names(row_df) == "marker"] <- marker_var

      for (cov_nm in cov_names) {
        row_df[[cov_nm]] <- dataEvent[[cov_nm]][i]
      }
      return(row_df)
    }
    do.call(rbind, rows_i)
  }

  obs_rows <- .sim_parallel_lapply(seq_len(n_id), .sim_obs_rows_for_id)
  obs_rows <- Filter(Negate(is.null), obs_rows)
  if (!length(obs_rows)) {
    cli::cli_abort(c(
      x = "No longitudinal observations remained after applying the schedule and event-time censoring.",
      i = "Use a denser {.arg times_obs} grid, reduce {.arg obs_time_noise_sd}, or disable {.arg censor_longitudinal_after_event}."
    ))
  }
  dataLong <- do.call(rbind, obs_rows)
  rownames(dataLong) <- NULL

  # ---- Mean structure from model matrices and sampled random effects
  dataLong_scaled <- dataLong
  dataLong_scaled[[time_var]] <- dataLong_scaled[[time_var]] / time_scale_internal
  X_long <- .mm(fixed_template, dataLong_scaled)
  Z_id_long <- .sim_rhs_matrix(id_templates, dataLong_scaled)
  Z_mk_long <- .sim_rhs_matrix(mk_templates, dataLong_scaled)
  Z_idm_long <- .sim_rhs_matrix(idm_templates, dataLong_scaled)

  id_index <- match(as.character(dataLong[[id_var]]), as.character(dataEvent[[id_var]]))
  marker_index <- match(as.character(dataLong[[marker_var]]), marker_levels)

  mu_long <- as.numeric(X_long %*% beta_long_internal)
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
  attr(dataLong, "JoiNMe_family_by_row") <- family_by_row

  sigma_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "sigma", 1.0), numeric(1))
  nu_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "nu", 4.0), numeric(1))
  phi_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi", 2.0), numeric(1))
  alpha_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "alpha", 0.0), numeric(1))
  kappa_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "kappa", 10.0), numeric(1))
  tau_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "tau", 0.5), numeric(1))
  trials_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "trials", 10L), numeric(1))

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
    if (param_name == "kappa") kappa_vec <- exp(eta_param)
    if (param_name == "tau") tau_vec <- stats::plogis(eta_param)
  }

  # A fixed family quantile is the final statistical specification for its
  # marker. Apply it after distributional regression so mixed simulations
  # mirror Stan: fixed markers use their declared constants, while remaining
  # skew-Laplace markers retain their family-level or row-specific values.
  fixed_tau_row <- use_tau_fixed[marker_index] == 1L
  tau_vec[fixed_tau_row] <- tau_fixed[marker_index[fixed_tau_row]]

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
        kappa = kappa_vec[r],
        tau = tau_vec[r],
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
  kappa_vec <- kappa_vec[ord_long]
  tau_vec <- tau_vec[ord_long]
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
    kappa = kappa_vec,
    tau = tau_vec,
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
    kappa = kappa_vec,
    tau = tau_vec,
    trials = as.integer(round(trials_vec)),
    coef = dist_coef_effective,
    eta = dist_eta_effective,
    design_cols = dist_design_cols,
    family_defaults = family_params,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels)
  )

  tmax <- max(dataEvent[[event_time_var]])
  time_scale_generation <- as.numeric(time_scale_internal %||% tmax)

  beta_scaled <- beta_long
  if (length(idx_time_beta) > 0L) {
    beta_scaled[idx_time_beta] <- beta_scaled[idx_time_beta] * time_scale_generation
  }

  tau_u_eff <- re_id_effective$sd
  if (length(idx_time_uid) > 0L) {
    tau_u_eff[idx_time_uid] <- tau_u_eff[idx_time_uid] * time_scale_generation
  }

  tau_v_eff <- re_marker_effective$sd
  if (length(idx_time_vmk) > 0L) {
    tau_v_eff[idx_time_vmk] <- tau_v_eff[idx_time_vmk] * time_scale_generation
  }

  marker_id_row_scale_eff <- rep(1.0, K_idm)
  if (length(idx_time_idm) > 0L) {
    marker_id_row_scale_eff[idx_time_idm] <- time_scale_generation
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
  kappa_family_truth <- .sim_family_param_truth("kappa", 10.0)
  # The fitted model declares a family-level tau only when at least one
  # skew-Laplace marker has not supplied its own fixed quantile.
  tau_family_truth <- if (any(skew_laplace_marker & use_tau_fixed == 0L)) {
    .sim_family_param_truth("tau", 0.5)
  } else {
    numeric(0)
  }
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
  marker_to_kappa_family <- setNames(match(family_names, names(kappa_family_truth), nomatch = 0L), marker_levels)
  marker_to_tau_family <- setNames(match(family_names, names(tau_family_truth), nomatch = 0L), marker_levels)
  marker_to_tau_family[use_tau_fixed == 1L] <- 0L

  stan_fit_truth <- list(
    beta = beta_long,
    beta_scaled = beta_scaled,
    time_scale_generation = as.numeric(time_scale_generation),
    time_scale_observed_max = as.numeric(tmax),
    idx_time_beta = as.integer(idx_time_beta),
    idx_time_uid = as.integer(idx_time_uid),
    idx_time_vmk = as.integer(idx_time_vmk),
    idx_time_idm = as.integer(idx_time_idm),
    gamma_w = beta_event,
    bs_gamma_c = baseline_hazard_truth$coefficients,
    shrinkage = shrinkage,
    fixed_marker_weights = as.integer(isTRUE(fixed_marker_weights)),
    marker_weights_base = marker_weights_base,
    z_marker_weights = as.numeric(t(z_marker_weight_sets)),
    marker_weights_eff = marker_weights_effective,
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
    alpha_cv_total = unname(assoc_coef_scalar[["cv_total"]] %||% NA_real_),
    alpha_cs_total = unname(assoc_coef_scalar[["cs_total"]] %||% NA_real_),
    alpha_cv_mean = unname(assoc_coef_scalar[["cv_mean"]] %||% NA_real_),
    alpha_cs_mean = unname(assoc_coef_scalar[["cs_mean"]] %||% NA_real_),
    alpha_cv_marker = unname(assoc_coef_scalar[["cv_marker"]] %||% NA_real_),
    alpha_cs_marker = unname(assoc_coef_scalar[["cs_marker"]] %||% NA_real_),
    alpha_corr = assoc_coef_corr,
    alpha_vcov = assoc_coef_vcov,
    sigma_family = sigma_family_truth,
    nu_family = nu_family_truth,
    phi_family = phi_family_truth,
    alpha_family = alpha_family_truth,
    kappa_family = kappa_family_truth,
    tau_family = tau_family_truth,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels),
    trials_family = trials_family_truth,
    cutpoints_ord = cutpoints_ord_truth,
    marker_to_sigma_family = marker_to_sigma_family,
    marker_to_nu_family = marker_to_nu_family,
    marker_to_phi_family = marker_to_phi_family,
    marker_to_alpha_family = marker_to_alpha_family,
    marker_to_kappa_family = marker_to_kappa_family,
    marker_to_tau_family = marker_to_tau_family
  )

  # Assemble a simulation-truth record in both the public statistical language
  # and the fitted Stan parameter language. The standardised draws are retained
  # because class recovery should be assessed at the level on which the mixture
  # is defined, not by classifying covariance-scaled effects after generation.
  mixture_truth <- mixture_configuration
  if (!is.null(mixture_truth)) {
    mixture_truth$formulaClass <-
      .mixture_specification$formulaClass # original membership formula request
    mixture_truth$standardised_draws <- list(
      subject = z_id,
      marker = z_marker
    ) # latent values to which component locations and scales may have been applied
    covariance_class_type <- intersect(
      c("corr", "vcov"),
      mixture_truth$class_type
    ) # selected public name for the covariance-regression latent matrix
    if (length(covariance_class_type) == 1L) {
      mixture_truth$standardised_draws[[covariance_class_type]] <-
        z_cov[, seq_len(M_cov), drop = FALSE]
    }
    stan_fit_truth$mix_probability <-
      unname(mixture_truth$probability)
    stan_fit_truth$mix_location <-
      unname(mixture_truth$location)
    stan_fit_truth$mix_scale <-
      unname(mixture_truth$scale)
    stan_fit_truth$mix_class_coefficient_subject <-
      unname(mixture_truth$coefficient$subject)
    stan_fit_truth$mix_class_coefficient_marker <-
      unname(mixture_truth$coefficient$marker)
    stan_fit_truth$posterior_class_subject <-
      unname(mixture_truth$allocation$subject)
    stan_fit_truth$posterior_class_marker <-
      unname(mixture_truth$allocation$marker)
  }

  family_truth <- list(
    by_marker = data.frame(
      marker = marker_levels,
      family = family_names,
      family_code = family_codes,
      marker_to_sigma_family = unname(marker_to_sigma_family),
      marker_to_nu_family = unname(marker_to_nu_family),
      marker_to_phi_family = unname(marker_to_phi_family),
      marker_to_alpha_family = unname(marker_to_alpha_family),
      marker_to_kappa_family = unname(marker_to_kappa_family),
      marker_to_tau_family = unname(marker_to_tau_family),
      use_tau_fixed = use_tau_fixed,
      tau_fixed = tau_fixed,
      stringsAsFactors = FALSE
    ),
    shared = stan_fit_truth[c(
      "sigma_family",
      "nu_family",
      "phi_family",
      "alpha_family",
      "kappa_family",
      "tau_family",
      "trials_family",
      "cutpoints_ord"
    )],
    defaults = family_params
  )

  true_params <- list(
    beta_long = beta_long,
    beta_event = beta_event,
    alpha_cv_total = if ("cv_total" %in% names(assoc_coef_vec)) assoc_coef_vec[["cv_total"]] else NA_real_,
    alpha_cs_total = if ("cs_total" %in% names(assoc_coef_vec)) assoc_coef_vec[["cs_total"]] else NA_real_,
    alpha_cv_mean = if ("cv_mean" %in% names(assoc_coef_vec)) assoc_coef_vec[["cv_mean"]] else NA_real_,
    alpha_cs_mean = if ("cs_mean" %in% names(assoc_coef_vec)) assoc_coef_vec[["cs_mean"]] else NA_real_,
    alpha_corr = assoc_coef_corr,
    alpha_vcov = assoc_coef_vcov,
    assoc = assoc,
    assoc_coefs = assoc_coef_vec,
    # `marker_weights` remains the convenient public truth alias and denotes
    # the effective weights that actually generated the event process.
    marker_weights = marker_weights_effective,
    marker_weights_raw = marker_weights_effective,
    marker_weights_base = marker_weights_base,
    marker_weights_latent = marker_weights_latent,
    marker_weights_eff = marker_weights_eff,
    marker_weights_by_term = marker_weights_eff,
    marker_weights_base_by_term = marker_weights_base_by_term,
    marker_weights_latent_by_term = marker_weights_latent_by_term,
    z_marker_weight_sets = z_marker_weight_sets,
    fixed_marker_weights = isTRUE(fixed_marker_weights),
    shrinkage = shrinkage,
    shrinkage_distribution = c(
      `0` = "student_t(6, 0, 1)",
      `1` = "double_exponential(0, 1)",
      `2` = "normal(0, 1)"
    )[[as.character(shrinkage)]],
    shared_marker_weights = isTRUE(marker_weight_spec$shared_marker_weights),
    link_names = link_names,
    quadrature_nodes = as.integer(gk_spec$n_gk),
    gk_nodes = gk_spec$nodes,
    gk_weights = gk_spec$weights,
    gk_rule = gk_spec$rule,
    transforms = transforms,
    basehaz = h0_fn,
    baseline_hazard = baseline_hazard_truth,
    formulaVCov = formulaVCov,
    formulaBasehaz = formulaBasehaz,
    beta_basehaz = baseline_hazard_truth$coefficients,
    family = family_truth,
    distributional_params = distributional_truth,
    stan_fit = stan_fit_truth,
    dist_coefs = dist_coefs,
    dist_re_params = re_params$dist %||% list(),
    dist_re = dist_re_effective %||% list(),
    re_params = re_params,
    re_structure = list(
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
      id = re_id_public,
      marker = re_marker_public,
      id_marker_cov = re_idm_public,
      id_marker_cov_latent = re_idm,
      id_marker_cov_scaled = re_idm_scaled
    ),
    re_draws_likelihood = list(
      id = re_id,
      marker = re_marker,
      id_marker_cov = re_idm_scaled
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
    mixture = mixture_truth,
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaDist = dist_formulas
  )

  marker_info <- list(
    names = marker_levels,
    families = families,
    family_codes = family_codes,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels)
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

  # ---- Step 5c: Build event rows automatically from the simulated risk-process
  # structure. Counting-process rows are emitted whenever split layout is
  # required by the data-generating mechanism.
  dataEvent_public <- dataEvent
  if (event_layout == "split") {
    interval_rows <- lapply(seq_len(n_id), function(i) {
      stop_i <- as.numeric(dataEvent[[event_time_var]][i])
      status_i <- as.integer(dataEvent[[event_var]][i])
      start_i <- as.numeric(entry_time_by_id[i])

      tv_internal <- numeric(0)
      if (length(tv_cov_specs) > 0L) {
        tv_internal <- unique(unlist(lapply(tv_cov_specs, function(spec_all) {
          breaks_i <- as.numeric(spec_all[[i]]$breaks %||% numeric(0))
          breaks_i[breaks_i > start_i & breaks_i < stop_i]
        }), use.names = FALSE))
      }
      internal <- times_obs[times_obs > start_i & times_obs < stop_i]
      breaks <- c(start_i, internal, tv_internal, stop_i)
      breaks <- sort(unique(as.numeric(breaks)))
      if (length(breaks) < 2L) {
        breaks <- c(start_i, stop_i)
      }
      int_start <- breaks[-length(breaks)]
      int_stop <- breaks[-1]
      int_status <- rep.int(0L, length(int_start))
      if (status_i == 1L && length(int_status) > 0L) {
        int_status[length(int_status)] <- 1L
      }

      out <- dataEvent[rep(i, length(int_start)), , drop = FALSE]
      if (length(tv_cov_specs) > 0L) {
        for (nm in names(tv_cov_specs)) {
          out[[nm]] <- vapply(int_start, function(tt) {
            .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], tt)
          }, numeric(1))
        }
      }
      out$time_start <- as.numeric(int_start)
      # Ensure that time_stop is always strictly greater than time_start to avoid zero-length intervals
      out$time_stop <- pmax(as.numeric(int_stop), out$time_start + 1e-9)
      out[[event_time_var]] <- pmax(as.numeric(int_stop), out$time_start + 1e-9)
      out[[event_var]] <- as.integer(int_status)
      out
    })
    dataEvent_public <- do.call(rbind, interval_rows)
    rownames(dataEvent_public) <- NULL
  } else {
    if (length(tv_cov_specs) > 0L) {
      for (nm in names(tv_cov_specs)) {
        dataEvent_public[[nm]] <- vapply(seq_len(nrow(dataEvent_public)), function(i) {
          .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], 0)
        }, numeric(1))
      }
    }
    dataEvent_public$time_start <- as.numeric(entry_time_by_id)
    # Ensure that time_stop is always strictly greater than time_start to avoid zero-length intervals
    dataEvent_public$time_stop <- 
      pmax(
        as.numeric(dataEvent_public[[event_time_var]]), 
        dataEvent_public$time_start + 1e-9)
  }

  list(
    dataLong = dataLong,
    dataEvent = dataEvent_public,
    truth = true_params,
    true_params = true_params,
    marker_info = marker_info,
    helpers = helpers,
    tmax = tmax
  )
}
