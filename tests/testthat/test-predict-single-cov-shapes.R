test_that("posterior_predict handles single-covariate hazard/corr shapes", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  set.seed(1421)
  sim <- simulate_joinme(
    n_id = 8,
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x2,
    formulaVCov = ~ x1,
    families = c("gaussian", "student_t", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 3, length.out = 4),
    assoc = c("cv_total", "corr"),
    seed = 1421
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x2,
    formulaVCov = ~ x1,
    families = c("gaussian", "student_t", "student_t"),
    dataLong = sim$dataLong,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total", "corr"),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 20,
      iter_sampling = 20,
      refresh = 0,
      seed = 1421
    )
  )

  last_time <- sim$dataLong |>
    tidytable::summarize(time_start = max(time), .by = id)

  ndE <- tidytable::left_join(sim$dataEvent, last_time, by = "id")

  pred <- posterior_predict(
    fit,
    newdataLong = tidytable::filter(sim$dataLong, id %in% c(1, 2)),
    newdataEvent = tidytable::filter(ndE, id %in% c(1, 2)),
    time_start = "time_start",
    time_horizon = 4,
    control = list(
      n_samples = 10,
      n_times = 10,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 20,
      iter_sampling = 20,
      refresh = 0
    )
  )

  expect_s3_class(pred, "JoinMeDynPred")
  expect_true(is.data.frame(pred$results$survival))
})
