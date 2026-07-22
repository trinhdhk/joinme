devtools::load_all(path = ".", quiet = TRUE)

compute_bfmi <- function(energy_mat) {
  apply(energy_mat, 2, function(e) {
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  })
}

run_fit <- function(sim, adapt_delta, seed_fit) {
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
      iter_warmup = 200,
      iter_sampling = 300,
      refresh = 100,
      adapt_delta = adapt_delta,
      max_treedepth = 14,
      seed = seed_fit,
      force_recompile = FALSE
    )
  )

  d <- fit$fit$diagnostic_summary()
  e <- fit$fit$sampler_diagnostics()[, , "energy__"]
  b <- compute_bfmi(e)

  data.frame(
    adapt_delta = adapt_delta,
    ebfmi_chain1 = d$ebfmi[1],
    ebfmi_chain2 = d$ebfmi[2],
    div_chain1 = d$num_divergent[1],
    div_chain2 = d$num_divergent[2],
    td_chain1 = d$num_max_treedepth[1],
    td_chain2 = d$num_max_treedepth[2],
    bfmi_manual_chain1 = b[1],
    bfmi_manual_chain2 = b[2]
  )
}

sim300 <- simulate_joinme(
  n_id = 300,
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  families = rep("gaussian", 3),
  n_obs_per_marker_per_id = 8,
  times_obs = seq(0, 6, length.out = 10),
  beta_long = c("(Intercept)" = 0.8, "time" = 0.6, "x1" = 0.5),
  beta_event = c("x2" = 0.1),
  assoc = c("cv_mean", "cv_marker"),
  assoc_coefs = c(cv_mean = 0.3, cv_marker = 0.1),
  marker_weights = c(1, 2, 1),
  fixed_marker_weights = TRUE,
  baseline_hazard = list(type = "weibull", shape = 1.1, scale = 10),
  seed = 20260303,
  time_cens = 10,
  use_mirai = TRUE,
  n_workers = 4
)

cat("Event rate:", mean(sim300$dataEvent$event), "\n")

res_090 <- run_fit(sim300, adapt_delta = 0.90, seed_fit = 20260310)
res_095 <- run_fit(sim300, adapt_delta = 0.95, seed_fit = 20260311)
res_099 <- run_fit(sim300, adapt_delta = 0.99, seed_fit = 20260312)

cat("\n===== TUNING COMPARISON =====\n")
print(rbind(res_090, res_095, res_099))
