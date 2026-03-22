suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, compile = FALSE)
})

set.seed(4101)

sim <- simulate_joinme(
  n_id = 3,
  families = rep("gaussian", 2),
  n_obs_per_marker_per_id = 3,
  times_obs = seq(0, 2, length.out = 4),
  quadrature_nodes = 7,
  seed = 4101,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.15)
)

formula_long <- y ~ 1 + splines::ns(time, df = 3) + x1 +
  (0 + splines::bs(time, df = 3) | id) +
  (0 + x1 + (0 + splines::ns(time, df = 3) | id) | marker)
formula_event <- survival::Surv(time, event) ~ 1 + x1 + x2

fit <- joinme(
  formulaLong = formula_long,
  dataLong = sim$dataLong,
  formulaEvent = formula_event,
  dataEvent = sim$dataEvent,
  formulaVCov = ~ 1,
  assoc = c("cv_total"),
  families = rep("gaussian", 2),
  transforms = list(cv_total = list(type = "identity")),
  control = list(
    engine = "cmdstanr",
    force_recompile = FALSE,
    quadrature_nodes = 7,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 10,
    iter_sampling = 10,
    seed = 4101,
    refresh = 0
  )
)

pred <- posterior_epred(
  fit,
  newdataLong = sim$dataLong,
  newdataEvent = sim$dataEvent,
  time_start = min(sim$dataLong$time),
  times = seq(0, max(sim$dataLong$time) + 0.5, length.out = 12),
  control = list(
    engine = "cmdstanr",
    force_recompile = FALSE,
    quadrature_nodes = 7,
    n_samples = 10,
    chains = 1,
    iter_warmup = 5,
    iter_sampling = 5,
    refresh = 0
  ),
  seed = 4101
)

cat("fit_class=", class(fit)[1], "\n", sep = "")
cat("pred_class=", class(pred)[1], "\n", sep = "")
cat("idx_time_beta_length=", length(fit$stan_data$idx_time_beta), "\n", sep = "")
cat("idx_time_uid_length=", length(fit$stan_data$idx_time_uid), "\n", sep = "")
cat("idx_time_idm_length=", length(fit$stan_data$idx_time_idm), "\n", sep = "")
pred_longitudinal <- pred$predictions$longitudinal
if (is.null(pred_longitudinal)) {
  pred_longitudinal <- pred$predictions$epred
}
cat("pred_longitudinal_rows=", nrow(pred_longitudinal), "\n", sep = "")