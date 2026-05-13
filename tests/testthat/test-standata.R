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

  expect_warning(
    sd <- joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      assoc = c("corr", "vcov")
    ),
    "corr.*ignored"
  )
  expect_equal(sd$assoc_corr, 0L)
  expect_equal(sd$assoc_vcov, 1L)
})

test_that("standata builds one marker-weight structure per active weighted association term", {
  set.seed(111)
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 3,
    n_t = 2,
    seed = 111,
    include_marker_only = TRUE
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  weights_by_term <- list(
    cv_total = c(m1 = 0.2, m2 = 0.5, m3 = 0.3),
    cs_total = c(m1 = -0.5, m2 = 0.1, m3 = 0.8),
    cv_marker = c(m1 = 1.0, m2 = 0.0, m3 = -1.0)
  )
  levels(sim$dataLong$marker) <- c("m1", "m2", "m3")

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total", "cs_total", "cv_marker"),
    marker_weights = weights_by_term,
    shared_marker_weights = FALSE
  )

  expect_equal(sd$shared_marker_weights, 0L)
  expect_equal(sd$n_marker_weight_sets, 3L)
  expect_equal(sd$marker_weight_set_cv_total, 1L)
  expect_equal(sd$marker_weight_set_cs_total, 2L)
  expect_equal(sd$marker_weight_set_cv_marker, 3L)
  expect_equal(sd$marker_weight_set_cs_marker, 0L)
  expect_equal(sd$marker_weights_by_term$cv_total, c(0.2, 0.5, 0.3), tolerance = 1e-8)
  expect_equal(sd$marker_weights_by_term$cs_total, c(-0.5, 0.1, 0.8), tolerance = 1e-8)
  expect_equal(sd$marker_weights_by_term$cv_marker, c(1.0, 0.0, -1.0), tolerance = 1e-8)
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

test_that("vcov dimension follows nested marker-by-id basis, not top-level id covariance", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 2,
    n_t = 2,
    seed = 206,
    include_marker_only = TRUE
  )

  sd_one <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + (1 | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("vcov"),
    transforms = joinme_tf(vcov = "identity")
  )

  expect_equal(sd_one$R_id, 2L)
  expect_equal(sd_one$Q_idm, 1L)
  expect_equal(sd_one$M_vcov_tf, 1L)

  sd_three <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("vcov"),
    transforms = joinme_tf(vcov = "identity")
  )

  expect_equal(sd_three$R_id, 2L)
  expect_equal(sd_three$Q_idm, 2L)
  expect_equal(sd_three$M_vcov_tf, 3L)
  expect_equal(sd_three$M_corr_tf, 1L)
})

test_that("standata marks top-level id double-bar terms as independent id covariance", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 3,
    n_t = 2,
    seed = 207,
    include_marker_only = TRUE
  )

  sd <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time || id) +
      (0 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("vcov"),
    transforms = joinme_tf(vcov = "identity")
  )

  expect_equal(sd$indep_id_re, 1L)
  expect_equal(sd$indep_marker_re, 0L)
  expect_equal(sd$indep_idmarker_cov, 0L)
  expect_equal(sd$allow_marker_crosscorr, 1L)
})

test_that("outer marker double-bar disables cross-correlation but keeps inner id covariance", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 3,
    n_t = 2,
    seed = 208,
    include_marker_only = TRUE
  )

  expect_warning(
    sd <- joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (x1 + (1 + time | id) || marker),
      dataLong = sim$dataLong,
      formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
      dataEvent = sim$dataEvent,
      assoc = c("corr"),
      allow_marker_crosscorr = 1L
    ),
    "cross-correlation.*set to 0"
  )

  expect_equal(sd$indep_marker_re, 1L)
  expect_equal(sd$indep_idmarker_cov, 0L)
  expect_equal(sd$allow_marker_crosscorr, 0L)
  expect_equal(sd$assoc_corr, 1L)
})

test_that("inner id double-bar rejects corr and keeps cross-correlation available", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 3,
    n_t = 2,
    seed = 209,
    include_marker_only = TRUE
  )

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (x1 + (1 + time || id) | marker),
      dataLong = sim$dataLong,
      formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
      dataEvent = sim$dataEvent,
      assoc = c("corr")
    ),
    "corr"
  )

  sd <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (x1 + (1 + time || id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("vcov")
  )

  expect_equal(sd$indep_marker_re, 0L)
  expect_equal(sd$indep_idmarker_cov, 1L)
  expect_equal(sd$allow_marker_crosscorr, 1L)
  expect_equal(sd$M_vcov_tf, sd$Q_idm)
})

test_that("functional vcov constants remain vector-shaped for Stan data", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 2,
    n_t = 2,
    seed = 205,
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
    transforms = joinme_tf(vcov = ~ log(1 + exp(x)))
  )

  expect_equal(sd$tf_mode_vcov, 1L)
  expect_equal(sd$n_const_vcov, 1L)
  expect_equal(as.numeric(sd$const_data_vcov), 1)

  sd_stan <- joinme:::.coerce_rstan_vectors(sd, c("const_data_vcov"))
  expect_true(is.array(sd_stan$const_data_vcov))
  expect_equal(dim(sd_stan$const_data_vcov), 1)
  expect_equal(as.numeric(sd_stan$const_data_vcov), 1)
})
