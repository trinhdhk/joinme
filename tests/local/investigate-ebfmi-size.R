devtools::load_all(path = ".", quiet = TRUE)

truth_cv_mean <- 0.3
truth_cv_marker <- 0.1

compute_bfmi <- function(energy_mat) {
  apply(energy_mat, 2, function(e) {
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  })
}

run_case <- function(n_id, seed_offset = 0L, adapt_delta = 0.75) {
  cat("\n==============================\n")
  cat("Running case n_id =", n_id, "\n")
  cat("==============================\n")

  sim <- simulate_joinme(
    n_id = n_id,
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x2,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 8,
    times_obs = seq(0, 6, length.out = 10),
    beta_long = c("(Intercept)" = 0.8, "time" = 0.6, "x1" = 0.5),
    beta_event = c("x2" = 0.1),
    assoc = c("cv_mean", "cv_marker"),
    assoc_coefs = c(cv_mean = truth_cv_mean, cv_marker = truth_cv_marker),
    marker_weights = c(1, 2, 1),
    fixed_marker_weights = TRUE,
    baseline_hazard = list(type = "weibull", shape = 1.1, scale = 10),
    seed = 20260225 + as.integer(seed_offset),
    time_cens = 10,
    use_mirai = TRUE,
    n_workers = 4
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_mean", "cv_marker"),
    families = rep("gaussian", 3),
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 6,
      iter_warmup = 300,
      iter_sampling = 500,
      refresh = 100,
      adapt_delta = adapt_delta,
      max_treedepth = 12,
      seed = 20260226 + as.integer(seed_offset),
      force_recompile = FALSE
    )
  )

  diag_summary <- fit$fit$diagnostic_summary()
  sampler_diag <- fit$fit$sampler_diagnostics()
  energy <- sampler_diag[, , "energy__"]
  bfmi_manual <- compute_bfmi(energy)

  summ <- summary(fit)$tables
  assoc_tab <- summ$assoc
  re_sd_tab <- summ$re_sd

  cat("n_event:", sum(sim$dataEvent$event), " / ", nrow(sim$dataEvent),
      " (rate=", round(mean(sim$dataEvent$event), 3), ")\n", sep = "")
  cat("Diagnostic summary:\n")
  print(diag_summary)

  cat("Manual BFMI by chain:\n")
  print(bfmi_manual)

  if (!is.null(assoc_tab)) {
    cat("Association estimates:\n")
    print(assoc_tab[assoc_tab$term %in% c("cv_mean", "cv_marker"), c("term", "Estimate", "Q2.5", "Q97.5")])
  }

  if (!is.null(re_sd_tab)) {
    cat("Random-effect SD summary (first 10 rows):\n")
    print(utils::head(re_sd_tab, 10))
  }

  invisible(list(
    n_id = n_id,
    fit = fit,
    diag_summary = diag_summary,
    bfmi_manual = bfmi_manual,
    n_event = sum(sim$dataEvent$event),
    event_rate = mean(sim$dataEvent$event)
  ))
}

res_200 <- run_case(200, seed_offset = 0L)
res_300 <- run_case(300, seed_offset = 1000L)

cat("\n\n===== BFMI COMPARISON =====\n")
print(data.frame(
  n_id = c(200, 300),
  ebfmi_chain1 = c(res_200$diag_summary$ebfmi[1], res_300$diag_summary$ebfmi[1]),
  ebfmi_chain2 = c(res_200$diag_summary$ebfmi[2], res_300$diag_summary$ebfmi[2]),
  bfmi_manual_chain1 = c(res_200$bfmi_manual[1], res_300$bfmi_manual[1]),
  bfmi_manual_chain2 = c(res_200$bfmi_manual[2], res_300$bfmi_manual[2]),
  event_rate = c(res_200$event_rate, res_300$event_rate)
))
