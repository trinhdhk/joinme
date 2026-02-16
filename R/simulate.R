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
  ifelse(x > 0, x - log1p(exp(-x)), log1p(exp(x)))
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
#' @param gamma_w Hazard covariates (Intercept, x1, x2).
#' @param alpha_cv_total Association coefficient for CV_total.
#' @param h0 Baseline hazard function h0(t).
#' @param t_admin Administrative censoring time.
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
  gamma_w = c("(Intercept)" = 0.0, "x1" = 0.3, "x2" = -0.2),
  alpha_cv_total = 0.6,
  h0 = weibull_h0(shape = 1.4, scale = 6.0),
  t_admin = 8.0,
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

  dist_formulas <- .normalize_formula_dist(formulaDist)

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

  # ---- marker-by-id latent z_w (id-specific)
  tau_w <- rep(1.0, Q_idm)
  Corr_w <- diag(Q_idm)
  if (Q_idm == 2) Corr_w <- matrix(c(1, 0.1, 0.1, 1), 2, 2)
  Lcorr_w <- t(chol(Corr_w))
  z_w_lat <- array(rnorm(n_id * D * Q_idm), dim = c(n_id, D, Q_idm))
  L_w <- diag(tau_w, Q_idm, Q_idm) %*% Lcorr_w
  z_w <- array(0.0, dim = c(n_id, D, Q_idm))
  for (i in seq_len(n_id)) {
    for (d in seq_len(D)) {
      z_w[i, d, ] <- as.numeric(L_w %*% z_w_lat[i, d, ])
    }
  }

  # ---- id-specific covariance L_i (diagonal)
  alpha_L <- rep(-0.2, Q_idm)
  lambda_L <- rep(0.3, Q_idm)
  tau_L <- 0.6
  z_L <- rnorm(n_id)
  u_L <- tau_L * z_L

  L_i <- array(0.0, dim = c(n_id, Q_idm, Q_idm))
  for (i in seq_len(n_id)) {
    for (q in seq_len(Q_idm)) {
      lp <- alpha_L[q] + lambda_L[q] * u_L[i]
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
  haz_row <- function(i) c("(Intercept)" = 1, x1 = dataEvent$x1[i], x2 = dataEvent$x2[i])

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
      return(list(time = t_admin, event = 0L, bracketing_failed = TRUE, U = U))
    }
    T <- uniroot(f, lower = br$lower, upper = br$upper)$root
    if (T > t_admin) {
      list(time = t_admin, event = 0L, bracketing_failed = FALSE, U = U)
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
      t_max_i <- min(dataEvent$time[i], t_admin)
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
        X_sigma <- .build_dist_matrix(dist_formulas$sigma, df)$X
        sig <- as.numeric(exp(X_sigma %*% beta_sigma))
      }
      if (length(sig) != length(mu)) {
        stop("sigma_fun must return a vector with length equal to times.")
      }
      nu_vec <- rep(nu, length(mu))
      if (!is.null(dist_formulas$nu) && !is.null(beta_nu)) {
        X_nu <- .build_dist_matrix(dist_formulas$nu, df)$X
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
    tau_w = tau_w,
    Lcorr_w = Lcorr_w,
    z_w_lat = z_w_lat,
    z_w = z_w,
    alpha_L = alpha_L,
    lambda_L = lambda_L,
    tau_L = tau_L,
    z_L = z_L,
    u_L = u_L,
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
#' Generalized simulator for `joinme` that mirrors the fitting syntax as closely as
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
#' @param formulaEvent Event/survival formula (same role as in `joinme()`).
#' @param formulaDist Optional distributional regression formulas (same role as in `joinme()`).
#' @param formulaAssoc Optional formula for association terms in hazard, e.g. `~ cv_total + vcov`.
#'   If provided, it overrides `assoc`.
#' @param n_id Number of subjects.
#' @param families Marker-specific family names.
#' @param marker_levels Optional marker names; defaults to `m1`, `m2`, ...
#' @param n_obs_per_marker_per_id Target number of observations per (id, marker).
#' @param times_obs Optional candidate observation time grid.
#' @param n_t Optional alias for `n_obs_per_marker_per_id` for compatibility.
#' @param seed RNG seed.
#' @param covariate_formulas Named or LHS formulas used to generate event-level covariates,
#'   e.g. `list(x1 ~ rnorm(n_id), x2 ~ rt(n_id, df = 5))`.
#' @param assoc Association components (same names as fit): `cv_total`, `cv_mean`,
#'   `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `vcov`.
#' @param assoc_coefs Named coefficients for association terms in the hazard.
#' @param beta_long Fixed-effect coefficients for `formulaLong` fixed part. If NULL,
#'   coefficients are randomly generated and named by model-matrix columns.
#' @param beta_event Survival baseline-covariate coefficients for `formulaEvent` RHS.
#'   If NULL, coefficients are randomly generated.
#' @param dist_coefs Named list of vectors for distributional regressions (`sigma`, `nu`,
#'   `phi`, `alpha`, `phi_beta`, `tau_sde`).
#' @param re_params Random-effects simulation controls.
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
#' @param t_admin Administrative censoring horizon.
#' @param eps_cs Finite-difference step for slope-type associations (`cs_*`).
#' @param integration_control Control list passed to `integrate()`.
#' @param root_control Root finding control for inverse-CDF sampling.
#' @param id_var,marker_var,time_var,y_var,event_time_var,event_var Column names aligned
#'   with `joinme_standata()` defaults.
#'
#' @return list(dataLong, dataEvent, true_params, marker_info, helpers, tmax)
#' @export
simulate_joinme <- function(
  formulaLong = y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
  formulaDist = NULL,
  formulaAssoc = NULL,
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
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6),
  beta_long = NULL,
  beta_event = NULL,
  dist_coefs = list(),
  re_params = list(
    id = list(sd = NULL, corr = NULL),
    marker = list(sd = NULL, corr = NULL),
    id_marker = list(sd = NULL, corr = NULL)
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
  t_admin = 8.0,
  eps_cs = 1e-3,
  integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1.0, t_max = 50.0, expand = 1.7, max_expand = 60L),
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
      out[names(user_coef)] <- as.numeric(user_coef)
    } else {
      out[seq_len(min(length(out), length(user_coef)))] <- as.numeric(user_coef)[seq_len(min(length(out), length(user_coef)))]
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

  .sim_draw_re_block <- function(n_group, K, cfg, label) {
    if (K <= 0) return(matrix(0.0, n_group, 0))
    sd_vec <- cfg$sd
    if (is.null(sd_vec)) sd_vec <- rep(0.5, K)
    if (length(sd_vec) == 1) sd_vec <- rep(sd_vec, K)
    if (length(sd_vec) != K) {
      cli::cli_abort(c(
        x = "Random-effect SD length mismatch for {label}.",
        i = "Expected length {K}, got {length(sd_vec)}."
      ))
    }
    corr <- cfg$corr
    if (is.null(corr)) corr <- diag(K)
    if (!is.matrix(corr) || any(dim(corr) != c(K, K))) {
      cli::cli_abort(c(
        x = "Random-effect correlation matrix mismatch for {label}.",
        i = "Expected a {K}x{K} matrix."
      ))
    }
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

  .sim_get_family_param <- function(fam_name, param_name, fallback) {
    val <- family_params[[fam_name]][[param_name]]
    if (is.null(val)) fallback else val
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
    unique(vars[vars %in% c("cv_total", "cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker", "vcov")])
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
      knots <- as.numeric(spec$knots %||% stats::quantile(seq(0, t_admin, length.out = 100), probs = c(0.25, 0.5, 0.75)))
      degree <- as.integer(spec$degree %||% 3L)
      coef <- spec$coef
      intercept <- as.logical(spec$intercept %||% TRUE)

      return(function(t) {
        t <- pmax(as.numeric(t), 0)
        if (basis_type == "ns") {
          B <- splines::ns(t, knots = knots, Boundary.knots = c(0, t_admin), intercept = intercept)
        } else {
          B <- splines::bs(t, knots = knots, Boundary.knots = c(0, t_admin), degree = degree, intercept = intercept)
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
  family_codes <- .parse_family(families)

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
  id_idx <- which(grp_names == id_var)
  id_rhs_list <- if (length(id_idx) > 0) .bar_terms_to_rhs_list(bars[id_idx]) else list()

  nested_terms <- .extract_nested_marker_terms(formulaLong, marker_var = marker_var, id_var = id_var)
  mk_rhs_list <- nested_terms$mk_rhs_list
  idm_rhs_list <- nested_terms$idm_rhs_list

  # ---- Build event-level frame and covariates
  dataEvent <- data.frame(id = seq_len(n_id), stringsAsFactors = FALSE)
  dataEvent <- .sim_eval_covariate_formulas(n_subjects = n_id, formulas = covariate_formulas, base_df = dataEvent)
  names(dataEvent)[names(dataEvent) == "id"] <- id_var

  # ---- Build prototype data for matrix column alignment
  prototype <- dataEvent[rep(1, D), , drop = FALSE]
  prototype[[marker_var]] <- factor(marker_levels, levels = marker_levels)
  prototype[[time_var]] <- rep(mean(range(times_obs)), D)

  X_proto <- .mm(fixed_rhs, prototype)
  Z_id_proto <- .sim_rhs_matrix(id_rhs_list, prototype)
  Z_mk_proto <- .sim_rhs_matrix(mk_rhs_list, prototype)
  Z_idm_proto <- .sim_rhs_matrix(idm_rhs_list, prototype)

  beta_long <- .sim_align_coef(colnames(X_proto), beta_long, sd_default = 0.35, intercept_default = 1.0)

  event_rhs <- stats::delete.response(stats::terms(formulaEvent))
  W_event <- stats::model.matrix(event_rhs, dataEvent)
  storage.mode(W_event) <- "double"
  beta_event <- .sim_align_coef(colnames(W_event), beta_event, sd_default = 0.25, intercept_default = 0.0)

  # ---- Draw random effects from user-configurable covariance structures
  K_id <- ncol(Z_id_proto)
  K_mk <- ncol(Z_mk_proto)
  K_idm <- ncol(Z_idm_proto)

  re_id <- .sim_draw_re_block(n_id, K_id, re_params$id %||% list(), "id")
  re_marker <- .sim_draw_re_block(D, K_mk, re_params$marker %||% list(), "marker")
  re_idm_flat <- .sim_draw_re_block(n_id * D, K_idm, re_params$id_marker %||% list(), "id_marker")

  re_idm <- array(0.0, dim = c(n_id, D, K_idm))
  if (K_idm > 0) {
    for (i in seq_len(n_id)) {
      for (d in seq_len(D)) {
        re_idm[i, d, ] <- re_idm_flat[(i - 1L) * D + d, ]
      }
    }
  }

  # ---- Association terms driving survival (syntax aligned with fit)
  assoc_from_formula <- .sim_assoc_from_formula(formulaAssoc)
  if (!is.null(assoc_from_formula)) assoc <- assoc_from_formula
  assoc <- unique(assoc)
  if (length(assoc) == 0) assoc <- "cv_total"

  assoc_coef_vec <- rep(0.0, length(assoc))
  names(assoc_coef_vec) <- assoc
  if (!is.null(names(assoc_coefs))) {
    common_assoc <- intersect(names(assoc_coefs), assoc)
    assoc_coef_vec[common_assoc] <- as.numeric(assoc_coefs[common_assoc])
  } else if (length(assoc_coefs) > 0) {
    assoc_coef_vec[seq_len(min(length(assoc), length(assoc_coefs)))] <- as.numeric(assoc_coefs)[seq_len(min(length(assoc), length(assoc_coefs)))]
  }

  # ---- Latent trajectory evaluators at arbitrary (id, marker, time)
  make_row_df <- function(i, d, t) {
    out <- dataEvent[i, , drop = FALSE]
    out[[time_var]] <- t
    out[[marker_var]] <- factor(marker_levels[d], levels = marker_levels)
    out
  }

  eta_components <- function(i, d, t) {
    row_df <- make_row_df(i, d, t)
    x_fix <- .mm(fixed_rhs, row_df)
    z_id <- .sim_rhs_matrix(id_rhs_list, row_df)
    z_mk <- .sim_rhs_matrix(mk_rhs_list, row_df)
    z_idm <- .sim_rhs_matrix(idm_rhs_list, row_df)

    fixed_part <- if (ncol(x_fix) > 0) as.numeric(x_fix %*% beta_long) else 0
    id_part <- if (ncol(z_id) > 0) sum(z_id[1, ] * re_id[i, ]) else 0
    mk_part <- if (ncol(z_mk) > 0) sum(z_mk[1, ] * re_marker[d, ]) else 0
    idm_part <- if (ncol(z_idm) > 0) sum(z_idm[1, ] * re_idm[i, d, ]) else 0

    list(
      fixed_id = fixed_part + id_part,
      marker_specific = mk_part + idm_part,
      total = fixed_part + id_part + mk_part + idm_part
    )
  }

  assoc_components <- function(i, t) {
    mu_total <- numeric(D)
    mu_mean <- numeric(D)
    mu_marker <- numeric(D)
    for (d in seq_len(D)) {
      part <- eta_components(i, d, t)
      mu_total[d] <- part$total
      mu_mean[d] <- part$fixed_id
      mu_marker[d] <- part$marker_specific
    }
    cv_total <- mean(mu_total)
    cv_mean <- mean(mu_mean)
    cv_marker <- mean(mu_marker)

    t_eps <- t + eps_cs
    mu_total_eps <- numeric(D)
    mu_mean_eps <- numeric(D)
    mu_marker_eps <- numeric(D)
    for (d in seq_len(D)) {
      part_eps <- eta_components(i, d, t_eps)
      mu_total_eps[d] <- part_eps$total
      mu_mean_eps[d] <- part_eps$fixed_id
      mu_marker_eps[d] <- part_eps$marker_specific
    }

    list(
      cv_total = cv_total,
      cv_mean = cv_mean,
      cv_marker = cv_marker,
      cs_total = (mean(mu_total_eps) - cv_total) / eps_cs,
      cs_mean = (mean(mu_mean_eps) - cv_mean) / eps_cs,
      cs_marker = (mean(mu_marker_eps) - cv_marker) / eps_cs,
      vcov = if (D > 1) stats::var(mu_total) else 0
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

  hazard_i <- function(i, t) {
    if (length(t) > 1) {
      return(vapply(t, function(tt) hazard_i(i, tt), numeric(1)))
    }
    comp <- assoc_components(i, t)
    assoc_lp <- 0.0
    for (nm in assoc) {
      assoc_lp <- assoc_lp + assoc_coef_vec[[nm]] * comp[[nm]]
    }
    h0_fn(t) * exp(eta_event_i[i] + assoc_lp)
  }

  cumhaz_i <- function(i, t) {
    if (t <= 0) return(0)
    out <- do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integration_control))
    as.numeric(out$value)
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
      return(list(time = t_admin, event = 0L, bracketing_failed = TRUE))
    }
    T_star <- stats::uniroot(f_root, lower = br$lower, upper = br$upper)$root
    if (T_star > t_admin) {
      list(time = t_admin, event = 0L, bracketing_failed = FALSE)
    } else {
      list(time = T_star, event = 1L, bracketing_failed = FALSE)
    }
  }

  # ---- Draw event/censoring times
  event_draws <- lapply(seq_len(n_id), draw_event_time)
  dataEvent[[event_time_var]] <- vapply(event_draws, `[[`, numeric(1), "time")
  dataEvent[[event_var]] <- vapply(event_draws, `[[`, integer(1), "event")

  # ---- Build longitudinal observation schedule conditional on event times
  obs_rows <- vector("list", n_id * D)
  idx_row <- 1L
  n_obs_target <- max(2L, as.integer(n_obs_per_marker_per_id))

  for (i in seq_len(n_id)) {
    obs_upper <- max(1e-8, min(dataEvent[[event_time_var]][i], t_admin))
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

      cov_names <- setdiff(colnames(dataEvent), c(id_var, event_time_var, event_var))
      for (cov_nm in cov_names) {
        row_df[[cov_nm]] <- dataEvent[[cov_nm]][i]
      }

      obs_rows[[idx_row]] <- row_df
      idx_row <- idx_row + 1L
    }
  }
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
    idm_effect <- matrix(0.0, nrow(dataLong), K_idm)
    for (r in seq_len(nrow(dataLong))) {
      idm_effect[r, ] <- re_idm[id_index[r], marker_index[r], ]
    }
    mu_long <- mu_long + rowSums(Z_idm_long * idm_effect)
  }

  # ---- Distributional parameters (family defaults + formulaDist overrides)
  family_by_row <- families[marker_index]

  sigma_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "sigma", 1.0), numeric(1))
  nu_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "nu", 4.0), numeric(1))
  phi_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi", 2.0), numeric(1))
  alpha_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "alpha", 0.0), numeric(1))
  phi_beta_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi_beta", 10.0), numeric(1))
  tau_sde_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "tau_sde", 0.5), numeric(1))
  trials_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "trials", 10L), numeric(1))

  dist_formulas <- .normalize_formula_dist(formulaDist)
  for (param_name in names(dist_formulas)) {
    X_param <- .build_dist_matrix(dist_formulas[[param_name]], dataLong)$X
    beta_param <- .sim_align_coef(colnames(X_param), dist_coefs[[param_name]], sd_default = 0.15, intercept_default = 0.0)
    eta_param <- as.numeric(X_param %*% beta_param)

    if (param_name == "sigma") sigma_vec <- exp(eta_param)
    if (param_name == "nu") nu_vec <- 2 + exp(eta_param)
    if (param_name == "phi") phi_vec <- exp(eta_param)
    if (param_name == "alpha") alpha_vec <- eta_param
    if (param_name == "phi_beta") phi_beta_vec <- exp(eta_param)
    if (param_name == "tau_sde") tau_sde_vec <- stats::plogis(eta_param)
  }

  # ---- Draw outcomes by marker-specific family
  y_out <- numeric(nrow(dataLong))
  for (r in seq_len(nrow(dataLong))) {
    fam_name <- family_by_row[r]
    fam_code <- .parse_family(fam_name)
    eta_r <- mu_long[r]

    if (fam_code == 11L) {
      cp <- family_params[[fam_name]]$cutpoints %||% c(-1, 1)
      y_out[r] <- .sim_sample_ordinal(eta = eta_r, cutpoints = cp)
    } else {
      y_out[r] <- .sample_from_family(
        n = 1,
        eta = eta_r,
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
  dataLong <- dataLong[order(dataLong[[id_var]], dataLong[[marker_var]], dataLong[[time_var]]), , drop = FALSE]
  rownames(dataLong) <- NULL

  true_params <- list(
    beta_long = beta_long,
    beta_event = beta_event,
    beta = beta_long,
    gamma_w = beta_event,
    alpha_cv_total = assoc_coef_vec[["cv_total"]] %||% NA_real_,
    assoc = assoc,
    assoc_coefs = assoc_coef_vec,
    baseline_hazard = baseline_hazard,
    formulaBasehaz = formulaBasehaz,
    beta_basehaz = beta_basehaz,
    dist_coefs = dist_coefs,
    re_params = re_params,
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
    vcov = function(i, t) assoc_components(i, t)$vcov,
    assoc_components = assoc_components,
    baseline_hazard = h0_fn,
    hazard = hazard_i,
    cumhaz = cumhaz_i,
    survival_prob = function(i, t) exp(-cumhaz_i(i, t))
  )

  tmax <- max(dataEvent[[event_time_var]])

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
