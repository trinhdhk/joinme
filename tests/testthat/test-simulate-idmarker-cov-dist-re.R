test_that("simulate_joinme supports id-specific id:marker covariance via formulaCorr", {
  set.seed(4401)
  sim <- simulate_joinme(
    n_id = 12,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4401,
    formulaCorr = ~ x1 + x2,
    re_params = list(
      id = list(sd = c(0.5, 0.2)),
      marker = list(sd = 0.3),
      id_marker_cov = list(
        latent = list(sd = c(0.5, 0.25)),
        alpha = c(-0.3, 0.15, -0.1),
        beta = matrix(0.1, nrow = 3, ncol = 2),
        lambda = 0.6,
        sd_u = 0.7,
        diag_link = "softplus"
      )
    )
  )

  raw1 <- sim$helpers$assoc_components_raw(1, 1.0)
  raw2 <- sim$helpers$assoc_components_raw(2, 1.0)

  expect_true(inherits(sim$truth$formulaCorr, "formula"))
  expect_equal(as.character(sim$truth$formulaCorr), as.character(~ x1 + x2))
  expect_true(length(raw1$corr_vals) > 0)
  expect_false(isTRUE(all.equal(raw1$corr_vals, raw2$corr_vals)))
})

test_that("simulate_joinme supports distributional random effects", {
  set.seed(4402)
  sim <- simulate_joinme(
    n_id = 10,
    families = c("student_t", "student_t"),
    n_obs_per_marker_per_id = 5,
    times_obs = seq(0, 2, length.out = 6),
    seed = 4402,
    formulaDist = list(sigma ~ 1 + (1 | id)),
    re_params = list(
      id = list(sd = c(0.4, 0.2)),
      marker = list(sd = 0.25),
      id_marker_cov = list(latent = list(sd = c(0.4, 0.2))),
      dist = list(sigma = list(sd = 0.8))
    )
  )

  expect_true(is.list(sim$truth$dist_re_params))
  expect_true("sigma" %in% names(sim$truth$dist_re_params))
  expect_equal(sim$truth$dist_re_params$sigma$sd, 0.8)
  expect_true(all(is.finite(sim$dataLong$y)))
})
