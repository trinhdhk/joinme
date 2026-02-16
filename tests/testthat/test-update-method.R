test_that("update.JoinMeFit refits with updated formulas", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")

  set.seed(404)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 4, length.out = 8),
    seed = 404,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  sim$dataLong$x2 <- rnorm(nrow(sim$dataLong))
  sim$dataEvent$x3 <- rnorm(nrow(sim$dataEvent))

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
    families = rep("student_t", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0,
      seed = 404
    )
  )

  fit2 <- update(
    fit,
    formulaLong = ~ . + x2,
    formulaEvent = ~ . + x3,
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0,
      seed = 405
    )
  )

  expect_s3_class(fit2, "JoinMeFit")
  expect_true("x2" %in% all.vars(fit2$formulaLong))
  expect_true("x3" %in% all.vars(fit2$formulaEvent))
})
