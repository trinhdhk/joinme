test_that("standata builds and RMS-normalizes marker weights", {
  set.seed(101)
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 3),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 3, length.out = 6),
    seed = 101,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  marker_weights <- c("m3" = 0.2, "m1" = 0.5, "m2" = -0.3)
  levels(sim$dataLong$marker) <- c("m1", "m2", "m3")

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total", "cs_marker"),
    marker_weights = marker_weights
  )

  expected_weights <- marker_weights[sd$marker_levels]
  expected_weights <- expected_weights / sqrt(mean(expected_weights^2))
  expect_equal(sqrt(mean(sd$marker_weights^2)), 1, tolerance = 1e-8)
  expect_equal(sd$marker_weights, as.numeric(expected_weights), tolerance = 1e-8)
  expect_equal(sd$indep_id_re, 0L)
  expect_equal(sd$indep_marker_re, 0L)
  expect_equal(sd$indep_idmarker_cov, 0L)
  expect_equal(sd$assoc_cv_total, 1L)
  expect_equal(sd$assoc_cs_marker, 1L)
  expect_equal(sd$assoc_cv_mean, 0L)
  expect_equal(sd$assoc_cv_marker, 0L)
  expect_equal(sd$assoc_cs_total, 0L)
  expect_equal(sd$assoc_cs_mean, 0L)
})

test_that("standata rejects all-zero marker weights", {
  set.seed(111)
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 3),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 3, length.out = 6),
    seed = 111,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  expect_error(
    joinme_standata(
      formulaLong = formulaLong,
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      assoc = c("cv_total"),
      marker_weights = c(0, 0, 0)
    ),
    "all zeros"
  )
})

test_that("standata resolves double-bar independence", {
  set.seed(303)
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 3, length.out = 6),
    seed = 303,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time || id) +
    (0 + x1 + (1 + time || id) || marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total")
  )

  expect_equal(sd$indep_id_re, 1L)
  expect_equal(sd$indep_marker_re, 1L)
  expect_equal(sd$indep_idmarker_cov, 1L)
})

test_that("standata validates formula requirements", {
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 3, length.out = 6),
    seed = 202,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + (1 + time | id),
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent
    ),
    "marker"
  )

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 | marker),
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      assoc = c("vcov")
    ),
    "vcov"
  )
})

test_that("standata supports vcov link and fixed tau", {
  set.seed(505)
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 3, length.out = 6),
    seed = 505,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = "vcov",
    vcov_diag_link = "exp"
  )

  expect_equal(sd$vcov_diag_link, 1L)

  sim_sde <- simulate_joinme(
    n_id = 3,
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 1, length.out = 5),
    seed = 506,
    families = "skew_double_exponential"
  )
  formulaLong_sde <- y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 | marker)

  sd_sde <- joinme_standata(
    formulaLong = formulaLong_sde,
    dataLong = sim_sde$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim_sde$dataEvent,
    families = "skew_double_exponential",
    tau_sde_fixed = 0.2
  )

  expect_equal(sd_sde$use_tau_sde_fixed, 1L)
  expect_equal(sd_sde$tau_sde_fixed, 0.2)

  sim_gauss <- simulate_joinme(
    n_id = 3,
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 1, length.out = 5),
    seed = 507,
    families = "gaussian"
  )
  formulaLong_gauss <- y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 | marker)

  expect_warning(
    joinme_standata(
      formulaLong = formulaLong_gauss,
      dataLong = sim_gauss$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim_gauss$dataEvent,
      families = "gaussian",
      tau_sde_fixed = 0.2
    ),
    "ignored"
  )
})
