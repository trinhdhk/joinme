test_that("standata builds marker weights", {
  set.seed(101)
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 3,
    n_t = 2,
    seed = 101,
    include_marker_only = TRUE
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  marker_weights <- c("m3" = 0.2, "m1" = 0.5, "m2" = 0.3)
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
  expect_equal(sd$marker_weights, as.numeric(expected_weights), tolerance = 1e-8)
  expect_equal(sd$assoc_cv_total, 1L)
  expect_equal(sd$assoc_cs_marker, 1L)
  expect_equal(sd$assoc_cv_mean, 0L)
  expect_equal(sd$assoc_cv_marker, 0L)
  expect_equal(sd$assoc_cs_total, 0L)
  expect_equal(sd$assoc_cs_mean, 0L)
})

test_that("standata validates formula requirements", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 2,
    n_t = 2,
    seed = 202,
    include_marker_only = TRUE
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
      assoc = c("corr")
    ),
    "corr"
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

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      assoc = c("corr", "vcov")
    ),
    "cannot be used together"
  )
})

test_that("standata builds vcov association metadata", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 2,
    n_t = 2,
    seed = 204,
    include_marker_only = TRUE
  )

  sd <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("vcov"),
    transforms = joinme_tf(vcov = "identity")
  )

  expect_equal(sd$assoc_corr, 0L)
  expect_equal(sd$assoc_vcov, 1L)
  expect_equal(sd$tf_mode_vcov, 0L)
})
