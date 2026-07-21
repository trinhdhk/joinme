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
  expect_equal(unname(sim$truth$alpha_corr), vc)
  expect_equal(unname(sim$truth$stan_fit$alpha_corr), vc)
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

test_that("simulate_joinme records scalar association truth with canonical output names", {
  expected <- c(
    alpha_cv_total = 0.1,
    alpha_cs_total = 0.3,
    alpha_cv_mean = 0.2,
    alpha_cs_mean = 0.4
  )
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 1, length.out = 3),
    seed = 441,
    assoc = c("cv_total", "cv_mean", "cs_total", "cs_mean"),
    assoc_coefs = c(cv_total = 0.1, cv_mean = 0.2, cs_total = 0.3, cs_mean = 0.4)
  )

  observed <- unlist(sim$truth$stan_fit[names(expected)], use.names = TRUE)
  expect_equal(observed, expected)
  expect_equal(unlist(sim$truth[names(expected)], use.names = TRUE), expected)
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
  expect_equal(unname(sim$truth$alpha_vcov), 2)
  expect_equal(unname(sim$truth$stan_fit$alpha_vcov), 2)
})

test_that("vcov helper returns lower-triangular Cholesky-factor entries", {
  Li <- matrix(c(
    0.7, 0.0,
    -0.2, 0.5
  ), nrow = 2, byrow = TRUE)

  sd2 <- sqrt(sum(Li[2, ]^2))
  expect_equal(joinme:::.assoc_corr_features_from_chol(Li), c(Li[2, 1] / sd2), tolerance = 1e-8)
  expect_equal(joinme:::.assoc_vcov_features_from_chol(Li), c(Li[2, 1] / sd2, 0.7, sd2), tolerance = 1e-8)
  expect_equal(joinme:::.assoc_vcov_features_from_chol(Li, diagonal_only = TRUE), c(0.7, sd2), tolerance = 1e-8)
})

test_that("simulate_joinme vcov uses correlation-factor and SD features", {
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
  expect_true(all(raw$vcov_vals[2:3] > 0))
  expect_true(abs(raw$vcov_vals[1]) <= 1)
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

test_that("simulate_joinme vcov uses effective time-scaled SD features", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    n_id = 4,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 3, length.out = 5),
    seed = 3316,
    assoc = "vcov",
    assoc_coefs = list(vcov = c(0.2, -0.1, 0.3)),
    transforms = list(vcov = list(type = "identity"))
  )

  raw <- sim$helpers$assoc_components_raw(1, 1.2)
  row_scale <- sim$truth$stan_fit$marker_id_row_scale_eff
  li_eff <- sweep(sim$truth$L_i[1, , ], 1L, row_scale, `*`)
  expected <- joinme:::.assoc_vcov_features_from_chol(li_eff)

  expect_equal(raw$vcov_vals, expected, tolerance = 1e-10)
  expect_gt(raw$vcov_vals[3], joinme:::.assoc_vcov_features_from_chol(sim$truth$L_i[1, , ])[3])
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

test_that("simulate_joinme warns and ignores corr when vcov is also requested", {
  expect_warning(
    sim <- simulate_joinme(
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
    "corr.*ignored"
  )

  expect_false(any(grepl("^corr\\[", names(sim$truth$assoc_coefs))))
  expect_true(any(grepl("^vcov\\[", names(sim$truth$assoc_coefs))))
})

test_that("simulate_joinme warns and diagonalizes id correlation under top-level || id", {
  expect_warning(
    sim <- simulate_joinme(
      formulaLong = y ~ 1 + time +
        (1 + time || id) +
        (0 + (1 | id) | marker),
      n_id = 6,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3317,
      re_params = list(
        id = list(
          sd = c(0.5, 0.3),
          corr = matrix(c(1, 0.4, 0.4, 1), nrow = 2)
        )
      )
    ),
    "Ignoring nonzero off-diagonal entries in `re_params\\$id\\$corr`"
  )

  expect_equal(sim$truth$stan_fit$Corr_u, diag(2), tolerance = 1e-10)
  expect_equal(sim$truth$stan_fit$Lcorr_u, diag(2), tolerance = 1e-10)
  expect_equal(sim$truth$stan_fit$Sigma_u, diag(c(0.5, 0.3)^2), tolerance = 1e-10)
})

test_that("simulate_joinme rejects corr when nested marker-by-id covariance is diagonal", {
  expect_error(
    simulate_joinme(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (x1 + (1 + time || id) | marker),
      formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
      n_id = 4,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3314,
      assoc = c("corr")
    ),
    "corr"
  )
})

test_that("simulate_joinme warns and truncates extra covariance-regression terms under nested || id", {
  expect_warning(
    sim <- simulate_joinme(
      formulaLong = y ~ 1 + time +
        (1 + time | id) +
        (0 + (1 + time || id) | marker),
      formulaEvent = survival::Surv(time, event) ~ 1,
      n_id = 6,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 3318,
      re_params = list(
        id_marker_cov = list(
          alpha = c(0.5, -0.2, -0.2),
          lambda = c(1, 0.5, 0.25)
        )
      )
    ),
    "re_params\\$id_marker_cov\\$(alpha|lambda)",
    all = TRUE
  )

  expect_equal(sim$truth$id_marker_cov_effective$alpha, c(0.5, -0.2), tolerance = 1e-10)
  expect_equal(sim$truth$id_marker_cov_effective$lambda, c(1, 0.5), tolerance = 1e-10)
})
