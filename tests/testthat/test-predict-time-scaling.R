test_that(".scale_draw_dependent_time_terms scales time-indexed draw components", {
  draws <- list(
    beta_fixed = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
    tau_id = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
    tau_marker = matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE)
  )

  stan_data <- list(
    idx_time_beta = c(2L),
    idx_time_uid = c(1L),
    idx_time_vmk = c(2L),
    idx_time_idm = c(1L)
  )

  out <- JoiNMe:::.scale_draw_dependent_time_terms(draws, stan_data, tmax = 10)

  expect_equal(out$beta_fixed[, 1], draws$beta_fixed[, 1])
  expect_equal(out$beta_fixed[, 2], draws$beta_fixed[, 2] * 10)
  expect_equal(out$beta_fixed[, 3], draws$beta_fixed[, 3])

  expect_equal(out$tau_id[, 1], draws$tau_id[, 1] * 10)
  expect_equal(out$tau_id[, 2], draws$tau_id[, 2])

  expect_equal(out$tau_marker[, 1], draws$tau_marker[, 1])
  expect_equal(out$tau_marker[, 2], draws$tau_marker[, 2] * 10)
})


test_that(".scale_draw_dependent_time_terms is no-op for tmax ~ 1", {
  draws <- list(
    beta_fixed = matrix(c(1, 2, 3), nrow = 1),
    tau_id = matrix(c(4, 5), nrow = 1),
    tau_marker = matrix(c(6, 7), nrow = 1)
  )

  stan_data <- list(
    idx_time_beta = c(1L, 3L),
    idx_time_uid = c(2L),
    idx_time_vmk = c(1L),
    idx_time_idm = c(2L)
  )

  out <- JoiNMe:::.scale_draw_dependent_time_terms(draws, stan_data, tmax = 1)

  expect_equal(out$beta_fixed, draws$beta_fixed)
  expect_equal(out$tau_id, draws$tau_id)
  expect_equal(out$tau_marker, draws$tau_marker)
})
