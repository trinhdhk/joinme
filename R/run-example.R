#' Example runner for JoiNMe
#'
#' @description
#' Runs a small end-to-end workflow: simulate data, build standata, fit the model,
#' and return both the fit and simulated data. This is intended as a quick smoke
#' test or reproducible demo with fixed settings.
#'
#' @return A list with `fit` and `sim` entries.
#' @export
run_joinme_example <- function() {
  # End-to-end smoke test: simulate -> fit -> summarise
  sim <- simulate_joinme(
    n_id = 50,
    families = rep("student_t", 10),
    n_obs_per_marker_per_id = 10,
    times_obs = seq(0, 8, length.out = 16),
    seed = 42
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)

  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaVCov = ~1,
    dataLong = sim$dataLong,
    dataEvent = sim$dataEvent,
    draws = 500,
    basehaz = "bs",
    n_knots = 5,
    basehaz_degree = 3,
    tau_spline = 0.4,
    eps_fd = 1e-3,
    assoc = c("cv_total"),
    transforms = joinme_tf(
      cv_total = "identity",
      corr = "identity"
    ),
    formulaDist = list(sigma = sigma ~ 1 + time),
    priors = joinme_priors(beta = list(scale = 2.5), alpha = list(scale = 1.0), lkj = 2),
    control = list(
      parallel_chains = 2,
      iter_warmup = 200,
      iter_sampling = 200,
      threads_per_chain = 8,
      seed = 123,
      init = 0.5,
      refresh = 50,
      adapt_delta = 0.75,
      max_treedepth = 12
    )
  )

  summary(fit, draws = 200)
  invisible(list(fit = fit, sim = sim))
}
