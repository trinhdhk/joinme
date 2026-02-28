devtools::load_all(path = ".", quiet = TRUE)

truth_cv_mean <- 0.30
truth_cs_mean <- 0.20

sim <- simulate_joinme(
  n_id = 300,
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  families = rep("gaussian", 5),
  n_obs_per_marker_per_id = 8,
  times_obs = seq(0, 6, length.out = 10),
  beta_long = c("(Intercept)" = 0.8, "time" = 0.6, "x1" = 0.5),
  beta_event = c("x2" = 0.15),
  assoc = c("cv_mean"),
  assoc_coefs = c(cv_mean = truth_cv_mean, cv_marker = truth_cs_mean),
  marker_weights = c(1, 0.7, 1.1, 0.9, 1.0),
  seed = 20260225,
  time_cens = 10,
  use_mirai = TRUE,
  n_workers = 4
)

fit <- joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  dataLong = sim$dataLong,
  formulaEvent = survival::Surv(time, event) ~ x2,
  dataEvent = sim$dataEvent,
  assoc = c("cv_mean"),
  families = rep("gaussian", 5),
  fixed_marker_weights = FALSE,
  control = list(
    engine = "cmdstanr",
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 6,
    iter_warmup = 200,
    iter_sampling = 300,
    refresh = 50,
    adapt_delta = 0.85,
    max_treedepth = 12,
    seed = 20260226
  )
)

tab <- summary(fit)$tables$assoc
sel <- tab[tab$term %in% c("cv_mean", "cv_marker"), c("term", "Estimate", "Q2.5", "Q97.5")]
truth <- data.frame(term = c("cv_mean", "cv_marker"), true = c(truth_cv_mean, truth_cs_mean))
out <- merge(truth, sel, by = "term", all.x = TRUE, sort = FALSE)
out$cover <- with(out, true >= Q2.5 & true <= Q97.5)
out$abs_err <- abs(out$Estimate - out$true)

print(data.frame(event_rate = mean(sim$dataEvent$event)))
print(out)

diag <- fit$fit$diagnostic_summary()
print(diag)
if (!is.null(diag$num_divergent) && sum(diag$num_divergent) > 0) {
  stop("Divergent transitions detected.")
}
