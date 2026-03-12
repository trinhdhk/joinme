test_that("all", {
  skip_on_cran()
  set.seed(42)
  sim <- simulate_joinme(
    n_id = 100,
    families = rep("student_t", 5),
    n_obs_per_marker_per_id = 10,
    times_obs = seq(0, 8, length.out = 16),
    seed = 42,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = "student_t",
    transforms = joinme_tf(cv_total = "identity"),
    fixed_marker_weights = FALSE,
    control = list(
      # engine = "rstan",
      parallel_chains = 2,
      iter_warmup = 100,
      iter_sampling = 100,
      refresh = 0,
      adapt_delta = 0.78,
      max_treedepth = 12,
      seed = 421
    )
  )

  expect_s3_class(fit, "JoinMeFit")
  sum_obj <- summary(fit)
  expect_true(!is.null(sum_obj$tables))

  pred <- rstantools::posterior_predict(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    time_start = max(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 1, length.out = 20),
    control = list(
      n_samples = 20,
      chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0
    )
  )

  expect_s3_class(pred, "JoinMeDynPred")
  p <- plot(pred, which = c("longitudinal", "survival"), combined = FALSE)
  expect_true(inherits(p, "ggplot") || is.list(p))
})
