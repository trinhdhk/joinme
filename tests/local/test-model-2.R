devtools::load_all(path = ".", quiet = TRUE)

truth_cv_mean <-  0.1
truth_cv_marker <- 0.3

sim <- simulate_joinme(
  n_id = 800,
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  families = rep("gaussian", 3),
  times_obs = seq(0, 6, length.out = 10),
  assoc = c("cv_mean", "cv_marker"),
  truth = jm_truth(
    longitudinal = c("(Intercept)" = 0.8, "time" = 0.6, "x1" = 0.5),
    survival = list(slope = c(x2 = 0.1)),
    assoc_coef = list(slope = c(cv_mean = truth_cv_mean, cv_marker = truth_cv_marker)),
    marker_weights = list(offset = c(-1, 2, 1), family = "constant"),
    basehaz = list(type = "weibull", shape = 1.1, scale = 10)
  ),
  transforms = joinme_tf(cv_marker = ~ expit(x)),
  seed = 202603,
  time_cens = 10,
  use_mirai = TRUE,
  n_workers = 5
)

run_fit <- function(label, alpha_prior, tau_spline) {
  cat("\n==== ", label, " ====\n", sep = "")

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_mean", "cv_marker"),
    families = rep("gaussian", 3),
    # priors = jm_prior(marker_weights = list(
    #   offset = sim$truth$marker_weights,
    #   family = "constant"
    # )),
    # priors = list(alpha = list(scale = alpha_prior)),
    # tau_spline = tau_spline,
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 6,
      iter_warmup = 1500,
      iter_sampling = 1000,
      refresh = 100,
      adapt_delta = 0.78,
      max_treedepth = 12,
      seed = 20260226,
      force_recompile = F
    )
  )

  tab <- summary(fit)$tables$assoc
  sel <- tab[tab$term %in% c("cv_mean", "cv_marker"), c("term", "Estimate", "Q2.5", "Q97.5")]
  truth <- data.frame(term = c("cv_mean", "cv_marker"), true = c(truth_cv_mean, truth_cv_marker))
  out <- merge(truth, sel, by = "term", all.x = TRUE, sort = FALSE)
  out$cover <- with(out, true >= Q2.5 & true <= Q97.5)
  out$abs_err <- abs(out$Estimate - out$true)

  print(data.frame(
    n_id = nrow(sim$dataEvent),
    n_event = sum(sim$dataEvent$event),
    event_rate = mean(sim$dataEvent$event)
  ))
  print(out)

  diag_summary <- fit$fit$diagnostic_summary()
  print(diag_summary)

  # ---- Diagnostics: compare key latent components by chain
  sd <- fit$stan_data
  get_draw_array <- function(vars) {
    posterior::as_draws_array(extract(
      fit,
      what = "raw",
      variable = vars,
      keep_chains = TRUE
    )$posterior_draws)
  }

  summ_by_chain <- function(arr) {
    # arr dims: iteration x chain x variable
    apply(arr, c(2, 3), mean)
  }

  i <- 1
  d <- 1

  beta_vars <- paste0("beta[", seq_len(sd$P), "]")
  u_vars <- paste0("u_id[", i, ",", seq_len(sd$R_id), "]")
  v_vars <- if (sd$R_mk > 0) paste0("v_marker[", d, ",", seq_len(sd$R_mk), "]") else character(0)
  w_vars <- if (sd$Q_idm > 0) paste0("w_idm[", i, ",", d, ",", seq_len(sd$Q_idm), "]") else character(0)
  alpha_vars <- c("alpha_cv_marker")
  cov_vars <- c(
    paste0("alpha_L[", seq_len(sd$Q_idm), "]"),
    paste0("lambda_L[", seq_len(sd$Q_idm), "]")
  )

  draws_key <- get_draw_array(c(beta_vars, u_vars, v_vars, w_vars, alpha_vars, cov_vars))
  means_key <- summ_by_chain(draws_key)
  print(list(chain_means = means_key))

  # ---- Approximate hazard at event time for subject 1 using posterior draws
  if (sd$K_event >= 1) {
    X_event <- sd$X_event_now[i, ]
    Z_id_event <- sd$Z_id_event_now[i, ]
    Z_mk_event <- if (sd$R_mk > 0) sd$Z_mk_event_now[i, ] else numeric(0)
    Z_idm_event <- if (sd$Q_idm > 0) sd$Z_idm_event_now[i, ] else numeric(0)
    W_event <- if (sd$p_w > 0) sd$W[i, ] else numeric(0)
    bs_event <- sd$Bs_event_c[i, ]

    beta_draws <- posterior::as_draws_matrix(get_draw_array(beta_vars))
    u_draws <- posterior::as_draws_matrix(get_draw_array(u_vars))
    v_draws <- if (length(v_vars) > 0) posterior::as_draws_matrix(get_draw_array(v_vars)) else NULL
    w_draws <- if (length(w_vars) > 0) posterior::as_draws_matrix(get_draw_array(w_vars)) else NULL
    alpha_draws <- posterior::as_draws_matrix(get_draw_array("alpha_cv_marker"))
    bs_gamma_draws <- posterior::as_draws_matrix(get_draw_array(paste0("bs_gamma_c[1,", seq_len(sd$Kbs), "]")))
    gamma_w_draws <- if (sd$p_w > 0) posterior::as_draws_matrix(get_draw_array(paste0("gamma_w[1,", seq_len(sd$p_w), "]"))) else NULL
    sigma_draws <- NULL
    if (!is.null(sd$n_family_sigma) && sd$n_family_sigma > 0) {
      sigma_draws <- posterior::as_draws_matrix(get_draw_array(paste0("sigma_family[", seq_len(sd$n_family_sigma), "]")))
    } else if ("sigma_y" %in% fit$fit$metadata()$variables) {
      sigma_draws <- posterior::as_draws_matrix(get_draw_array("sigma_y"))
    }

    cvm_S <- as.numeric(beta_draws %*% X_event + u_draws %*% Z_id_event)
    cvk_S <- rep(0, nrow(beta_draws))
    if (sd$R_mk > 0 && !is.null(v_draws)) cvk_S <- cvk_S + as.numeric(v_draws %*% Z_mk_event)
    if (sd$Q_idm > 0 && !is.null(w_draws)) cvk_S <- cvk_S + as.numeric(w_draws %*% Z_idm_event)

    w_eff <- sd$marker_weights[d]
    cv_marker_S <- w_eff * cvk_S
    eta_assoc_S <- as.numeric(alpha_draws) * cv_marker_S

    log_h0_S <- as.numeric(bs_gamma_draws %*% bs_event)
    eta_w <- if (!is.null(gamma_w_draws)) as.numeric(gamma_w_draws %*% W_event) else 0

    log_haz_draws <- log_h0_S + eta_w + eta_assoc_S
    haz_draws <- exp(log_haz_draws)
    truth_time <- sim$dataEvent$time[i]
    truth_haz <- sim$helpers$hazard(i, truth_time)
    truth_cv_marker <- sim$helpers$cv_marker(i, truth_time)
    truth_log_h0 <- log(sim$helpers$baseline_hazard(truth_time))

    intercept_draws <- bs_gamma_draws[, 1]
    non_intercept_draws <- if (ncol(bs_gamma_draws) > 1) bs_gamma_draws[, -1, drop = FALSE] else NULL
    log_h0_intercept_only <- intercept_draws * bs_event[1]

    bs_event_col_sd <- apply(sd$Bs_event_c, 2, stats::sd)
    const_cols <- which(bs_event_col_sd < 1e-8 | !is.finite(bs_event_col_sd))
    const_vals <- if (length(const_cols) > 0) sd$Bs_event_c[1, const_cols] else numeric(0)

    surv_truth <- sim$helpers$survival_prob(i, truth_time)
    cumhaz_draws <- posterior::as_draws_matrix(get_draw_array(paste0("cumhaz_event[", i, "]")))
    surv_prob_draws <- exp(-cumhaz_draws)
    cumhaz_truth <- sim$helpers$cumhaz(i, truth_time)

    print(list(
      cv_marker_draw_mean = mean(cv_marker_S),
      cv_marker_truth = truth_cv_marker,
      log_hazard_draw_mean = mean(log_haz_draws),
      hazard_draw_mean = mean(haz_draws),
      hazard_truth = truth_haz,
      cumhaz_truth = cumhaz_truth,
      cumhaz_draw_mean = mean(cumhaz_draws),
      surv_prob_truth = surv_truth,
      surv_prob_draw_mean = mean(surv_prob_draws),
      surv_prob_draw_q025 = stats::quantile(surv_prob_draws, 0.025),
      surv_prob_draw_q975 = stats::quantile(surv_prob_draws, 0.975),
      log_h0_truth = truth_log_h0,
      sigma_draw_mean = if (!is.null(sigma_draws)) mean(sigma_draws) else NA_real_,
      log_h0_mean = mean(log_h0_S),
      log_h0_min = min(log_h0_S),
      log_h0_max = max(log_h0_S),
      log_h0_intercept_only_mean = mean(log_h0_intercept_only),
      eta_w_mean = mean(eta_w),
      alpha_cv_marker_mean = mean(alpha_draws),
      bs_gamma_mean = if (!is.null(bs_gamma_draws)) mean(bs_gamma_draws) else NA_real_,
      bs_gamma_intercept_mean = mean(intercept_draws),
      bs_gamma_intercept_sd = stats::sd(intercept_draws),
      bs_gamma_non_intercept_mean = if (!is.null(non_intercept_draws)) mean(non_intercept_draws) else NA_real_,
      bs_event_first = bs_event[1],
      bs_event_range = range(bs_event),
      bs_event_intercept_range = range(sd$Bs_event_c[, 1]),
      bs_gk_intercept_range = range(sd$Bs_gk_c[, , 1]),
      bs_event_constant_cols = const_cols,
      bs_event_constant_vals = const_vals
    ))
  }

  if (!is.null(diag_summary$num_divergent) && sum(diag_summary$num_divergent) > 0) {
    message("Divergent transitions detected in ", label, ".")
  }

  invisible(fit)
}

run_fit("option_1_relaxed_prior", alpha_prior = 5, tau_spline = 1.0)
run_fit("option_2_weaker_smoothing", alpha_prior = 5, tau_spline = 5.0)
