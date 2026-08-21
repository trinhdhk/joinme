test_that("simulate_joinme keeps event covariates and uses random beta_event defaults when omitted", {
  set.seed(5001)
  sim <- simulate_joinme(
    n_id = 30,
    families = rep("gaussian", 3),
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_mean"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_mean = 0.25))),
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
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_mean"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_mean = 0.25)),
      survival = list(slope = c(x1 = 0, x2 = 0))
    ),
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
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.25))),
    seed = 5003
  )

  expect_gt(stats::sd(sim$dataEvent$x1), 0)
  expect_gt(stats::sd(sim$dataEvent$x2), 0)
  expect_false(all(abs(sim$truth$beta_event) < 1e-12))
})

test_that("simulate_joinme keeps spline-expanded event term labels in beta_event truth", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 | id),
    formulaEvent = survival::Surv(time, event) ~ splines::bs(x1, df = 3) + x2,
    n_id = 4,
    families = rep("gaussian", 2),
    time_cens = 1,
    times_obs = seq(0, 1, length.out = 3),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0))),
    seed = 5004,
    use_mirai = FALSE
  )

  expect_equal(
    names(sim$truth$beta_event),
    c(
      "splines::bs(x1, df = 3)1",
      "splines::bs(x1, df = 3)2",
      "splines::bs(x1, df = 3)3",
      "x2"
    )
  )
})
