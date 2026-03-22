test_that("simulate_joinme corr uses off-diagonal correlation features", {
  vc <- c(0.4, -0.1, 0.2)
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time + I(time^2) | id) | marker),
    n_id = 5,
    families = c("gaussian", "gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 3, length.out = 5),
    seed = 3301,
    assoc = "corr",
    assoc_coefs = list(corr = vc),
    transforms = list(corr = list(type = "identity"))
  )

  raw <- sim$helpers$assoc_components_raw(1, 1.2)

  expect_length(raw$corr_vals, 3)
  expect_true(all(is.finite(raw$corr_vals)))
  expect_true(all(abs(raw$corr_vals) <= 1))
  expect_equal(raw$corr, sum(vc * raw$corr_vals), tolerance = 1e-10)

  vc_names <- paste0("corr[", seq_along(vc), "]")
  expect_equal(unname(sim$truth$assoc_coefs[vc_names]), vc)
})

test_that("simulate_joinme corr coefficient parsing truncates to M_corr", {
  expect_warning(
    sim <- simulate_joinme(
      n_id = 4,
      families = c("gaussian", "gaussian", "gaussian"),
      n_obs_per_marker_per_id = 4,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3302,
      assoc = "corr",
      assoc_coefs = c(corr = 0.7, corr2 = -0.5),
      transforms = list(corr = list(type = "identity"))
    ),
    "Ignoring extra `corr` association coefficients"
  )

  raw <- sim$helpers$assoc_components_raw(1, 0.9)

  expect_length(raw$corr_vals, 1)
  expect_equal(unname(sim$truth$assoc_coefs["corr[1]"]), 0.7)
  expect_equal(raw$corr, 0.7 * raw$corr_vals[1], tolerance = 1e-10)
})

test_that("simulate_joinme warns when extra vcov coefficients are supplied", {
  expect_warning(
    sim <- simulate_joinme(
      formulaLong = y ~ 1 + time +
        (1 + time | id) +
        (0 + (1 | id) | marker),
      n_id = 4,
      families = c("gaussian", "gaussian", "gaussian"),
      n_obs_per_marker_per_id = 4,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3303,
      assoc = "vcov",
      assoc_coefs = list(vcov = c(2, -1, -1)),
      transforms = list(vcov = list(type = "identity"))
    ),
    "Ignoring extra `vcov` association coefficients"
  )

  raw <- sim$helpers$assoc_components_raw(1, 0.9)

  expect_length(raw$vcov_vals, 1)
  expect_equal(unname(sim$truth$assoc_coefs["vcov[1]"]), 2)
})

test_that("vcov helper returns lower-triangular Cholesky-factor entries", {
  Li <- matrix(c(
    0.7, 0.0,
    -0.2, 0.5
  ), nrow = 2, byrow = TRUE)

  expect_equal(joinme:::.assoc_vcov_features_from_chol(Li), c(0.7, -0.2, 0.5))
  expect_equal(joinme:::.assoc_vcov_features_from_chol(Li, diagonal_only = TRUE), c(0.7, 0.5))
})

test_that("simulate_joinme vcov uses lower-triangular Cholesky-factor features", {
  vc <- c(0.2, -0.1, 0.3)
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    n_id = 5,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 3, length.out = 5),
    seed = 3311,
    assoc = "vcov",
    assoc_coefs = list(vcov = vc),
    transforms = list(vcov = list(type = "identity"))
  )

  raw <- sim$helpers$assoc_components_raw(1, 1.2)

  expect_length(raw$vcov_vals, 3)
  expect_true(all(is.finite(raw$vcov_vals)))
  expect_true(all(raw$vcov_vals[c(1, 3)] > 0))
  expect_equal(raw$vcov, sum(vc * raw$vcov_vals), tolerance = 1e-10)

  vc_names <- paste0("vcov[", seq_along(vc), "]")
  expect_equal(unname(sim$truth$assoc_coefs[vc_names]), vc)
})

test_that("simulate_joinme uses component-specific covariance latents", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    families = c("gaussian", "gaussian", "gaussian"),
    n_id = 40,
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 4),
    assoc = "vcov",
    assoc_coefs = list(vcov = c(2, -1, -1)),
    transforms = joinme_tf(vcov = ~ log(1 + exp(x))),
    re_params = list(
      id_marker_cov = list(
        alpha = c(-0.5, 0.5, 1),
        lambda = c(0.1, -0.15, 0.2)
      )
    ),
    seed = 2401,
    use_mirai = FALSE
  )

  vals <- t(vapply(seq_len(nrow(sim$dataEvent)), function(i) {
    sim$helpers$assoc_components_raw(i, 1.2)$vcov_vals_tf
  }, numeric(3)))
  corr_vals <- stats::cor(vals)

  expect_equal(dim(sim$truth$id_marker_cov_effective$z), c(nrow(sim$dataEvent), 3L))
  expect_lt(max(abs(corr_vals[upper.tri(corr_vals)])), 0.9999)
})

test_that("simulate_joinme zero-references corr functional transforms at raw zero", {
  cc <- c(0.4)
  sim <- simulate_joinme(
    n_id = 5,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 3, length.out = 5),
    seed = 3314,
    assoc = "corr",
    assoc_coefs = list(corr = cc),
    transforms = list(corr = list(type = "functional", expr = ~ log(1 + exp(x))))
  )

  raw <- sim$helpers$assoc_components_raw(1, 1.2)
  expected_tf <- log1p(exp(raw$corr_vals)) - log1p(exp(0))

  expect_equal(raw$corr_vals_tf, expected_tf, tolerance = 1e-10)
  expect_equal(raw$corr, sum(cc * expected_tf), tolerance = 1e-10)
})

test_that("simulate_joinme zero-references vcov functional transforms at raw zero", {
  vc <- c(0.2, -0.1, 0.3)
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    n_id = 5,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 3, length.out = 5),
    seed = 3315,
    assoc = "vcov",
    assoc_coefs = list(vcov = vc),
    transforms = list(vcov = list(type = "functional", expr = ~ log(1 + exp(x))))
  )

  raw <- sim$helpers$assoc_components_raw(1, 1.2)
  expected_tf <- log1p(exp(raw$vcov_vals)) - log1p(exp(0))

  expect_equal(raw$vcov_vals_tf, expected_tf, tolerance = 1e-10)
  expect_equal(raw$vcov, sum(vc * expected_tf), tolerance = 1e-10)
})

test_that("simulate_joinme vcov respects diagonal-only marker-by-id independence under mirai", {
  skip_if_not_installed("mirai")

  expect_warning(
    sim <- simulate_joinme(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time || id) | marker),
      n_id = 4,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3312,
      assoc = "vcov",
      assoc_coefs = list(vcov = c(0.4, 0.6, 0.9)),
      transforms = list(vcov = list(type = "identity")),
      n_workers = 2,
      use_mirai = TRUE
    ),
    "Ignoring extra `vcov` association coefficients"
  )

  raw <- sim$helpers$assoc_components_raw(1, 0.9)

  expect_length(raw$vcov_vals, 2)
  expect_true(all(raw$vcov_vals > 0))
  expect_equal(unname(sim$truth$assoc_coefs[c("vcov[1]", "vcov[2]")]), c(0.4, 0.6))
  expect_equal(raw$vcov, sum(c(0.4, 0.6) * raw$vcov_vals), tolerance = 1e-10)
})

test_that("simulate_joinme rejects corr and vcov together", {
  expect_error(
    simulate_joinme(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      n_id = 4,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3313,
      assoc = c("corr", "vcov")
    ),
    "cannot be used together"
  )
})
