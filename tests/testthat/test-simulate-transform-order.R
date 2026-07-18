test_that("simulate_joinme applies cv_total transform before marker-weight averaging", {
  eps_cs <- 1e-3
  sim <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1201,
    assoc = c("cv_total", "cs_total"),
    assoc_coefs = c(cv_total = 0.0, cs_total = 0.0),
    marker_weights = c(2, 1),
    fixed_marker_weights = TRUE,
    transforms = list(
      cv_total = list(type = "functional", expr = ~ x^2)
    ),
    eps_cs = eps_cs
  )

  t0 <- 0.9
  comp <- sim$helpers$assoc_components(1, t0)
  raw_now <- sim$helpers$assoc_components_raw(1, t0)
  raw_eps <- sim$helpers$assoc_components_raw(1, t0 + eps_cs)
  w <- sim$truth$marker_weights
  D <- length(w)

  expected_cv_total <- sum(w * raw_now$marker_values$mu_total^2) / D
  expected_cv_total_old <- (sum(w * raw_now$marker_values$mu_total) / D)^2
  expected_cs_total <- (sum(w * raw_eps$marker_values$mu_total) / D -
    sum(w * raw_now$marker_values$mu_total) / D) / eps_cs

  expect_equal(comp$cv_total, expected_cv_total, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(comp$cv_total, expected_cv_total_old, tolerance = 1e-8)))
  expect_equal(comp$cs_total, expected_cs_total, tolerance = 1e-6)
})

test_that("simulate_joinme applies cv_marker transform before marker-weight averaging", {
  sim <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1202,
    assoc = c("cv_marker"),
    assoc_coefs = c(cv_marker = 0.0),
    marker_weights = c(1.5, 0.5, 2.0),
    fixed_marker_weights = TRUE,
    transforms = list(
      cv_marker = list(type = "functional", expr = ~ softplus(x))
    )
  )

  t0 <- 1.1
  comp <- sim$helpers$assoc_components(2, t0)
  raw_now <- sim$helpers$assoc_components_raw(2, t0)
  w <- sim$truth$marker_weights
  D <- length(w)

  expected_cv_marker <- sum(w * log1p(exp(raw_now$marker_values$mu_marker))) / D
  expected_cv_marker_old <- log1p(exp(sum(w * raw_now$marker_values$mu_marker) / D))

  expect_equal(comp$cv_marker, expected_cv_marker, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(comp$cv_marker, expected_cv_marker_old, tolerance = 1e-8)))
})

test_that("simulate_joinme applies cs_total transform before marker-weight averaging", {
  eps_cs <- 1e-3
  sim <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1203,
    assoc = c("cs_total"),
    assoc_coefs = c(cs_total = 0.0),
    marker_weights = c(2, -1),
    fixed_marker_weights = TRUE,
    transforms = list(
      cs_total = list(type = "functional", expr = ~ x^2)
    ),
    eps_cs = eps_cs
  )

  t0 <- 0.8
  comp <- sim$helpers$assoc_components(1, t0)
  raw_now <- sim$helpers$assoc_components_raw(1, t0)
  raw_eps <- sim$helpers$assoc_components_raw(1, t0 + eps_cs)
  w <- sim$truth$marker_weights
  D <- length(w)

  cs_total_marker <- (raw_eps$marker_values$mu_total - raw_now$marker_values$mu_total) / eps_cs
  expected_cs_total <- sum(w * (cs_total_marker^2)) / D
  expected_cs_total_old <- (sum(w * cs_total_marker) / D)^2

  expect_equal(comp$cs_total, expected_cs_total, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(comp$cs_total, expected_cs_total_old, tolerance = 1e-8)))
})

test_that("simulate_joinme applies cs_marker transform before marker-weight averaging", {
  eps_cs <- 1e-3
  sim <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1204,
    assoc = c("cs_marker"),
    assoc_coefs = c(cs_marker = 0.0),
    marker_weights = c(1.5, 0.5, 2.0),
    fixed_marker_weights = TRUE,
    transforms = list(
      cs_marker = list(type = "functional", expr = ~ softplus(x))
    ),
    eps_cs = eps_cs
  )

  t0 <- 1.0
  comp <- sim$helpers$assoc_components(2, t0)
  raw_now <- sim$helpers$assoc_components_raw(2, t0)
  raw_eps <- sim$helpers$assoc_components_raw(2, t0 + eps_cs)
  w <- sim$truth$marker_weights
  D <- length(w)

  cs_marker_vals <- (raw_eps$marker_values$mu_marker - raw_now$marker_values$mu_marker) / eps_cs
  expected_cs_marker <- sum(w * log1p(exp(cs_marker_vals))) / D
  expected_cs_marker_old <- log1p(exp(sum(w * cs_marker_vals) / D))

  expect_equal(comp$cs_marker, expected_cs_marker, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(comp$cs_marker, expected_cs_marker_old, tolerance = 1e-8)))
})

test_that("simulate_joinme rejects expit spline knots outside [0, 1]", {
  expect_error(
    simulate_joinme(
      n_id = 3,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 1205,
      assoc = c("cv_total"),
      assoc_coefs = c(cv_total = 0),
      transforms = list(
        cv_total = list(
          type = "ispline_expit",
          knots = c(-0.1, 0.5, 1.1),
          coeff = c(0, 0.5, 1)
        )
      )
    ),
    "must lie on the expit scale"
  )
})
