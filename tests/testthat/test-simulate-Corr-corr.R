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
  expect_equal(raw$corr, sum(abs(vc) * raw$corr_vals), tolerance = 1e-10)

  vc_names <- paste0("corr[", seq_along(vc), "]")
  expect_equal(unname(sim$truth$assoc_coefs[vc_names]), abs(vc))
})

test_that("simulate_joinme corr coefficient parsing truncates to M_corr", {
  sim <- simulate_joinme(
    n_id = 4,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 4),
    seed = 3302,
    assoc = "corr",
    assoc_coefs = c(corr = 0.7, corr2 = -0.5),
    transforms = list(corr = list(type = "identity"))
  )

  raw <- sim$helpers$assoc_components_raw(1, 0.9)

  expect_length(raw$corr_vals, 1)
  expect_equal(unname(sim$truth$assoc_coefs["corr[1]"]), 0.7)
  expect_equal(raw$corr, 0.7 * raw$corr_vals[1], tolerance = 1e-10)
})

