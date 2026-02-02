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

#' Stable softplus
#' @keywords internal
.softplus <- function(x) {
  ifelse(x > 20, x, log1p(exp(x)))
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

#' Simulate joint model data with mixed families
#'
#' @description
#' Simulate multivariate longitudinal + survival data with marker-specific families.
#'
#' @param n_id Number of subjects.
#' @param n_obs_per_marker_per_id Observations per marker per subject.
#' @param families Character vector of family names (one per marker).
#' @param times_obs Observation times (e.g., seq(0, 5, length.out = 10)).
#' @param seed Random seed.
#' @param family_params List of additional parameters per family.
#'
#' @return list(dataLong, dataEvent, true_params, marker_info)
#' @export
simulate_joinme <- function(
  n_id = 50,
  n_obs_per_marker_per_id = 8,
  families = c("gaussian", "student_t", "binomial"),
  times_obs = seq(0, 5, length.out = 8),
  seed = 42,
  family_params = list(
    gaussian = list(sigma = 1.0),
    student_t = list(sigma = 1.5, nu = 4),
    binomial = list(trials = 10),
    poisson = list(),
    negbin2 = list(phi = 2),
    skew_normal = list(sigma = 1.0, alpha = 0),
    double_exponential = list(sigma = 1.0),
    skew_double_exponential = list(sigma = 1.0, tau_sde = 0.5),
    beta = list(phi_beta = 10),
    cumulative_logit = list(cutpoints = c(-1, 1))
  )
) {
  set.seed(seed)

  D <- length(families)
  family_codes <- integer(D)
  for (d in seq_len(D)) {
    family_codes[d] <- .parse_family(families[d])
  }

  marker_names <- if (D <= 3) {
    c("m1", "m2", "m3")[seq_len(D)]
  } else {
    paste0("marker_", seq_len(D))
  }

  beta <- c(0.5, -0.2, 0.1)
  sigma_b <- c(0.8, 0.4)
  tau_w <- 0.6

  obs_list <- vector("list", D * n_id)
  idx_obs <- 1

  for (i in seq_len(n_id)) {
    b_i <- c(rnorm(1, 0, sigma_b[1]), rnorm(1, 0, sigma_b[2]))

    for (d in seq_len(D)) {
      w_id <- rnorm(1, 0, tau_w)
      n_obs <- sample(n_obs_per_marker_per_id - 2, 1) + 2
      times <- sort(sample(times_obs, n_obs, replace = FALSE))

      for (t in times) {
        x1 <- rnorm(1, 0, 1)
        x2 <- rnorm(1, 0, 1)
        eta <- beta[1] + b_i[1] + (beta[2] + b_i[2]) * t + w_id + 0.3 * x1 - 0.2 * x2

        fam <- families[d]
        if (fam == "gaussian") {
          sig <- family_params$gaussian$sigma %||% 1.0
          y <- rnorm(1, eta, sig)
        } else if (fam == "student_t") {
          sig <- family_params$student_t$sigma %||% 1.5
          nu <- family_params$student_t$nu %||% 4
          y <- stats::rt(1, df = nu) * sig + eta
        } else if (fam == "binomial") {
          trials <- family_params$binomial$trials %||% 10
          p <- stats::plogis(eta)
          y <- stats::rbinom(1, size = trials, prob = p)
        } else if (fam == "bernoulli") {
          p <- stats::plogis(eta)
          y <- stats::rbinom(1, size = 1, prob = p)
        } else if (fam == "poisson") {
          mu <- exp(eta)
          y <- stats::rpois(1, lambda = mu)
        } else if (fam == "negbin2") {
          mu <- exp(eta)
          phi <- family_params$negbin2$phi %||% 2
          y <- stats::rnbinom(1, size = phi, mu = mu)
        } else if (fam == "skew_normal") {
          sig <- family_params$skew_normal$sigma %||% 1.0
          y <- rnorm(1, eta, sig)
        } else if (fam == "double_exponential") {
          sig <- family_params$double_exponential$sigma %||% 1.0
          sign <- sample(c(-1, 1), size = 1)
          y <- eta + sign * rexp(1, rate = 1 / sig)
        } else if (fam == "skew_double_exponential") {
          sig <- family_params$skew_double_exponential$sigma %||% 1.0
          tau <- family_params$skew_double_exponential$tau_sde %||% 0.5
          u <- runif(1)
          e <- rexp(1, rate = 1)
          y <- if (u < tau) eta + e * sig / tau else eta - e * sig / (1 - tau)
        } else if (fam == "beta") {
          phi_beta <- family_params$beta$phi_beta %||% 10
          mu <- plogis(eta)
          shape1 <- pmax(mu * phi_beta, 1e-6)
          shape2 <- pmax((1 - mu) * phi_beta, 1e-6)
          y <- rbeta(1, shape1 = shape1, shape2 = shape2)
        } else if (fam == "cumulative_logit") {
          # Simple 3-category ordinal outcome using fixed cutpoints
          cutpoints <- family_params$cumulative_logit$cutpoints %||% c(-1, 1)
          p1 <- plogis(cutpoints[1] - eta)
          p2 <- plogis(cutpoints[2] - eta) - p1
          p3 <- 1 - plogis(cutpoints[2] - eta)
          y <- sample(1:3, size = 1, prob = c(p1, p2, p3))
        } else {
          cli::cli_abort(c(x = "Unknown family: {fam}.", i = "Check the families argument."))
        }

        obs_list[[idx_obs]] <- data.frame(
          id = i,
          marker = marker_names[d],
          time = t,
          y = y,
          x1 = x1,
          x2 = x2,
          stringsAsFactors = FALSE
        )
        idx_obs <- idx_obs + 1
      }
    }
  }

  dataLong <- do.call(rbind, obs_list)

  dataEvent <- data.frame(
    id = seq_len(n_id),
    time = runif(n_id, min = min(times_obs), max = max(times_obs) + 2),
    event = stats::rbinom(n_id, size = 1, prob = 0.7),
    x1 = rnorm(n_id),
    x2 = rnorm(n_id),
    stringsAsFactors = FALSE
  )

  list(
    dataLong = dataLong,
    dataEvent = dataEvent,
    true_params = list(beta = beta, sigma_b = sigma_b, tau_w = tau_w),
    marker_info = list(names = marker_names, families = families, family_codes = family_codes)
  )
}
