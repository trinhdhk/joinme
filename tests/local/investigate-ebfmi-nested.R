devtools::load_all(path = ".", quiet = TRUE)

compute_bfmi <- function(energy_mat) {
  apply(energy_mat, 2, function(e) {
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  })
}

fit_case <- function(dataLong, dataEvent, label, seed_fit) {
  cat("\n------------------------------\n")
  cat("Fitting case:", label, "\n")
  cat("------------------------------\n")

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = dataLong,
    formulaEvent = survival::Surv(time, event) ~ x2,
    dataEvent = dataEvent,
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
      adapt_delta = 0.75,
      max_treedepth = 12,
      seed = seed_fit,
      force_recompile = FALSE
    )
  )

  diag_summary <- fit$fit$diagnostic_summary()
  energy <- fit$fit$sampler_diagnostics()[, , "energy__"]
  bfmi_manual <- compute_bfmi(energy)
  assoc_tab <- summary(fit)$tables$assoc

  out <- list(
    label = label,
    n_id = nrow(dataEvent),
    event_rate = mean(dataEvent$event),
    n_event = sum(dataEvent$event),
    ebfmi = diag_summary$ebfmi,
    n_div = diag_summary$num_divergent,
    n_td = diag_summary$num_max_treedepth,
    bfmi_manual = bfmi_manual,
    assoc = assoc_tab[assoc_tab$term %in% c("cv_mean", "cv_marker"), c("term", "Estimate", "Q2.5", "Q97.5")]
  )

  print(list(
    n_id = out$n_id,
    n_event = out$n_event,
    event_rate = out$event_rate,
    ebfmi = out$ebfmi,
    n_div = out$n_div,
    n_td = out$n_td,
    bfmi_manual = out$bfmi_manual
  ))
  print(out$assoc)

  invisible(out)
}

set.seed(20260302)

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
  baseline_hazard = list(type = "weibull", shape = 1.1, scale = 10),
  seed = 20260303,
  time_cens = 10,
  use_mirai = TRUE,
  n_workers = 4
)

ids200 <- sort(unique(sim300$dataEvent$id))[1:200]
dataEvent200 <- sim300$dataEvent[sim300$dataEvent$id %in% ids200, , drop = FALSE]
dataLong200 <- sim300$dataLong[sim300$dataLong$id %in% ids200, , drop = FALSE]

res200 <- fit_case(dataLong200, dataEvent200, "nested-200", seed_fit = 20260304)
res300 <- fit_case(sim300$dataLong, sim300$dataEvent, "nested-300", seed_fit = 20260305)

cat("\n===== NESTED COMPARISON =====\n")
print(data.frame(
  label = c(res200$label, res300$label),
  n_id = c(res200$n_id, res300$n_id),
  event_rate = c(res200$event_rate, res300$event_rate),
  ebfmi_chain1 = c(res200$ebfmi[1], res300$ebfmi[1]),
  ebfmi_chain2 = c(res200$ebfmi[2], res300$ebfmi[2]),
  n_div_chain1 = c(res200$n_div[1], res300$n_div[1]),
  n_div_chain2 = c(res200$n_div[2], res300$n_div[2]),
  bfmi_manual_chain1 = c(res200$bfmi_manual[1], res300$bfmi_manual[1]),
  bfmi_manual_chain2 = c(res200$bfmi_manual[2], res300$bfmi_manual[2])
))
