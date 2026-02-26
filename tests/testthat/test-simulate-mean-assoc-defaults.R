test_that("simulate_joinme keeps event covariates and uses random beta_event defaults when omitted", {
  set.seed(5001)
  sim <- simulate_joinme(
    n_id = 30,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0.25),
    seed = 5001
  )

  expect_gt(stats::sd(sim$dataEvent$x1), 0)
  expect_gt(stats::sd(sim$dataEvent$x2), 0)
  expect_false(all(abs(sim$truth$beta_event) < 1e-12))
})

test_that("simulate_joinme does not freeze event covariates when beta_event is provided", {
  set.seed(5002)
  sim <- simulate_joinme(
    n_id = 30,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0.25),
    beta_event = c("(Intercept)" = 0.1, "x1" = 0, "x2" = 0),
    seed = 5002
  )

  expect_gt(stats::sd(sim$dataEvent$x1), 0)
  expect_gt(stats::sd(sim$dataEvent$x2), 0)
})

test_that("simulate_joinme keeps event covariates for non-mean-only assoc", {
  set.seed(5003)
  sim <- simulate_joinme(
    n_id = 30,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.25),
    seed = 5003
  )

  expect_gt(stats::sd(sim$dataEvent$x1), 0)
  expect_gt(stats::sd(sim$dataEvent$x2), 0)
  expect_false(all(abs(sim$truth$beta_event) < 1e-12))
})
