test_that("standata maps weighted id grouping to subject and RE weights", {
  set.seed(2401)
  sim <- simulate_joinme(
    n_id = 5,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 2401,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  id_w <- seq_len(nrow(sim$dataEvent)) / 2 + 0.5
  sim$dataEvent$id_w <- id_w
  sim$dataLong$id_w <- id_w[match(sim$dataLong$id, sim$dataEvent$id)]
  sim$dataLong$marker_w <- ifelse(sim$dataLong$marker == levels(sim$dataLong$marker)[1], 1, 2)

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | weighted(id, weights = id_w)) +
    (0 + x1 + (1 + time | weighted(id, weights = id_w)) | weighted(marker, weights = marker_w))
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total")
  )

  expect_equal(sd$subject_weights, id_w, tolerance = 1e-8)
  expect_equal(sd$re_weight_id, id_w, tolerance = 1e-8)
  expect_equal(sd$re_weight_idm, id_w, tolerance = 1e-8)
  expect_equal(sd$re_weight_L, id_w, tolerance = 1e-8)
  expect_equal(sd$re_weight_marker, c(1, 2), tolerance = 1e-8)
})

test_that("nested id weighting controls corr latent id weights", {
  set.seed(2403)
  sim <- simulate_joinme(
    n_id = 5,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 2403,
    assoc = c("cv_total", "corr"),
    assoc_coefs = c(cv_total = 0.2, corr = 0.1)
  )

  id_w_nested <- seq_len(nrow(sim$dataEvent)) / 3 + 0.7
  sim$dataEvent$id_w_nested <- id_w_nested
  sim$dataLong$id_w_nested <- id_w_nested[match(sim$dataLong$id, sim$dataEvent$id)]

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | weighted(id, weights = id_w_nested)) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total", "corr")
  )

  expect_equal(sd$re_weight_id, rep(1, nrow(sim$dataEvent)), tolerance = 1e-8)
  expect_equal(sd$re_weight_idm, id_w_nested, tolerance = 1e-8)
  expect_equal(sd$re_weight_L, id_w_nested, tolerance = 1e-8)
  expect_equal(sd$subject_weights, id_w_nested, tolerance = 1e-8)
})

test_that("standata builds weighted distributional random-effects weights", {
  set.seed(2402)
  sim <- simulate_joinme(
    n_id = 4,
    families = c("gaussian", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 2402,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  id_w <- seq_len(nrow(sim$dataEvent)) + 1
  sim$dataEvent$id_w <- id_w
  sim$dataLong$id_w <- id_w[match(sim$dataLong$id, sim$dataEvent$id)]

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    families = c("gaussian", "student_t"),
    formulaDist = list(
      sigma ~ 1 + (1 | weighted(id, weights = id_w))
    ),
    assoc = c("cv_total")
  )

  expect_equal(sd$n_re_sigma, 1L)
  expect_equal(as.numeric(sd$re_weight_sigma[[1]][1:sd$G_sigma[1]]), id_w, tolerance = 1e-8)
})

test_that("weighted grouping requires weights= argument (weight= rejected)", {
  set.seed(2404)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 2404,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  sim$dataEvent$id_w <- runif(nrow(sim$dataEvent), 0.8, 1.2)
  sim$dataLong$id_w <- sim$dataEvent$id_w[match(sim$dataLong$id, sim$dataEvent$id)]

  formulaLong_bad <- y ~ 1 + time + x1 +
    (1 + time | weighted(id, weight = id_w)) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  expect_error(
    joinme_standata(
      formulaLong = formulaLong_bad,
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      assoc = c("cv_total")
    ),
    "weights"
  )
})

