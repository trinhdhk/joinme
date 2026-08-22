devtools::load_all(path = ".", quiet = TRUE)

compute_bfmi <- function(energy_mat) {
  apply(energy_mat, 2, function(e) {
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  })
}

fit_one <- function(sim, marker_weight_family, label, seed_fit) {
  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_mean", "cv_marker"),
    families = rep("gaussian", 3),
    priors = jm_prior(marker_weights = list(
      offset = sim$truth$marker_weights,
      family = marker_weight_family
    )),
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 6,
      iter_warmup = 200,
      iter_sampling = 300,
      refresh = 100,
      adapt_delta = 0.9,
      max_treedepth = 14,
      seed = seed_fit,
      force_recompile = FALSE
    )
  )

  d <- fit$fit$diagnostic_summary()
  e <- fit$fit$sampler_diagnostics()[, , "energy__"]
  b <- compute_bfmi(e)
  assoc_tab <- summary(fit)$tables$assoc

  out <- data.frame(
    label = label,
    marker_weight_family = marker_weight_family,
    ebfmi_chain1 = d$ebfmi[1],
    ebfmi_chain2 = d$ebfmi[2],
    div_chain1 = d$num_divergent[1],
    div_chain2 = d$num_divergent[2],
    bfmi_manual_chain1 = b[1],
    bfmi_manual_chain2 = b[2]
  )

  print(out)
  print(assoc_tab[assoc_tab$term %in% c("cv_mean", "cv_marker"), c("term", "Estimate", "Q2.5", "Q97.5")])
  out
}

sim <- simulate_joinme(
  n_id = 300,
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  families = rep("gaussian", 3),
  times_obs = seq(0, 6, length.out = 10),
  assoc = c("cv_mean", "cv_marker"),
  truth = jm_truth(
    longitudinal = c("(Intercept)" = 0.8, "time" = 0.6, "x1" = 0.5),
    survival = list(slope = c(x2 = 0.1)),
    assoc_coef = list(slope = c(cv_mean = 0.3, cv_marker = 0.1)),
    marker_weights = list(offset = c(1, 2, 1), family = "constant"),
    basehaz = list(type = "weibull", shape = 1.1, scale = 10)
  ),
  seed = 20260303,
  time_cens = 10,
  use_mirai = TRUE,
  n_workers = 4
)

cat("event_rate:", mean(sim$dataEvent$event), "\n")

res_free <- fit_one(sim, marker_weight_family = "normal", label = "fitted_weights", seed_fit = 20260320)
res_fixed <- fit_one(sim, marker_weight_family = "constant", label = "constant_weights", seed_fit = 20260321)

cat("\n===== WEIGHT IDENTIFICATION COMPARISON =====\n")
print(rbind(res_free, res_fixed))
