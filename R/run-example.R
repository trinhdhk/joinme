#' Example runner for joinme
#'
#' @description
#' Runs an simulation + standata build + fit.
#' @param ... Actually no. Currently everything is fixed. But will support input later :p
#'
#' @export
run_joinme_example <- function() {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 50,
    D = 10,
    n_t = 10,
    seed = 42,
    include_marker_only = TRUE,
    formulaDist = list(sigma = sigma ~ 1 + time),
    beta_sigma = c("(Intercept)" = log(0.4), "time" = 0.15)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)

  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaVcov = ~1,
    dataLong = sim$dataLong,
    dataEvent = sim$dataEvent,
    draws = 500,
    basehaz = "bs",
    n_knots = 5,
    basehaz_degree = 3,
    tau_spline = 0.4,
    eps_fd = 1e-3,
    assoc = c("cv_total"),
    transforms = list(
      cv_total = list(type = "identity"),
      vcov = list(type = "identity")
    ),
    formulaDist = list(sigma = sigma ~ 1 + time),
    priors = list(beta = list(scale = 2.5), alpha = list(scale = 1.0), lkj = 2),
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
