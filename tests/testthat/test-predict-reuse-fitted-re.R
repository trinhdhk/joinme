test_that("fitted subject identifiers use the standata ordering", {
  object <- structure(list(
    stan_data = list(n_id = 3L),
    call = quote(joinme(id_var = "participant")),
    dataLong = data.frame(participant = c("b", "a", "c", "a"))
  ), class = "JoiNMeFit")

  expect_equal(
    .fitted_subject_indices(object, c("c", "a"), "participant"),
    c(c = 3L, a = 1L)
  )
  expect_error(
    .fitted_subject_indices(object, "new", "participant"),
    "Unknown participant value"
  )
})

test_that("fitted subject identifiers retain numeric ordering", {
  object <- structure(list(
    stan_data = list(n_id = 3L),
    call = quote(joinme(id_var = "participant")),
    dataLong = data.frame(participant = c(10, 2, 1))
  ), class = "JoiNMeFit")

  expect_identical(
    unname(.fitted_subject_indices(object, c(1, 2, 10), "participant")),
    1:3
  )
})

test_that("posterior prediction wrappers expose fitted-effect reuse", {
  expect_true("reuse_fitted_re" %in% names(formals(predict.JoiNMeFit)))
  expect_true("reuse_fitted_re" %in% names(formals(posterior_linpred.JoiNMeFit)))
  expect_true("reuse_fitted_re" %in% names(formals(posterior_epred.JoiNMeFit)))
  expect_true("reuse_fitted_re" %in% names(formals(posterior_predict.JoiNMeFit)))
  expect_true("reuse_fitted_re" %in% names(formals(predict.JoiNMeMixFit)))
})

test_that("realised fitted random effects remain paired by posterior draw", {
  fitted_draws <- cbind(
    `u_id[1,1]` = c(1, 2),
    `u_id[2,1]` = c(3, 4),
    `v_marker[1,1]` = c(5, 6),
    `v_marker[2,1]` = c(7, 8),
    `z_w[1,1,1]` = c(0.1, 0.2),
    `z_w[1,2,1]` = c(0.3, 0.4),
    `z_w[2,1,1]` = c(0.5, 0.6),
    `z_w[2,2,1]` = c(0.7, 0.8),
    `L_i[1,1,1]` = c(1.1, 1.2),
    `L_i[2,1,1]` = c(1.3, 1.4)
  ) # transformed posterior quantities stored by the fitting programme
  draw_list <- list(
    beta_fixed = matrix(c(0, 0), ncol = 1L),
    fitted_random_effect_draws = fitted_draws
  )
  stan_data <- list(
    D = 2L,
    R_id = 1L,
    R_mk = 1L,
    Q_idm = 1L,
    indep_idmarker_cov = 1L
  )

  reused <- .fitted_random_effect_stan_data(
    draws_list = draw_list,
    stan_data = stan_data,
    fitted_subject_index = 2L,
    reuse_fitted_re = TRUE
  )
  expect_identical(reused$reuse_fitted_re, 1L)
  expect_equal(as.numeric(reused$fitted_u_id), c(3, 4))
  expect_equal(as.numeric(reused$fitted_v_marker[, 1, 1]), c(5, 6))
  expect_equal(as.numeric(reused$fitted_z_w[, 2, 1]), c(0.7, 0.8))
  expect_equal(as.numeric(reused$fitted_L_i[, 1, 1]), c(1.3, 1.4))
  expect_true(all(reused$z_u == 0))

  neutral <- .fitted_random_effect_stan_data(
    draws_list = draw_list,
    stan_data = stan_data,
    reuse_fitted_re = FALSE
  )
  expect_identical(neutral$reuse_fitted_re, 0L)
  expect_true(all(neutral$fitted_u_id == 0))
  expect_equal(as.numeric(neutral$fitted_L_i[, 1, 1]), c(1, 1))
})

test_that("ordinary and mixture fitpred programmes share the fitted-effect contract", {
  required_names <- c(
    "reuse_fitted_re",
    "fitted_u_id",
    "fitted_v_marker",
    "fitted_z_w",
    "fitted_L_i",
    "z_u",
    "z_v",
    "z_w_lat",
    "z_L"
  )
  for (program in c("joinme_fitpred", "joinme_mix_fitpred")) {
    stan_file <- .get_stan_file(program)
    required_data <- .stan_data_names(stan_file)
    expect_true(file.exists(stan_file))
    expect_true(all(required_names %in% required_data))
  }
})

test_that("neutral covariance latents retain every packed lower-triangle coordinate", {
  draw_list <- list(
    beta_fixed = matrix(0, nrow = 2L, ncol = 1L),
    fitted_random_effect_draws = matrix(0, nrow = 2L, ncol = 0L)
  )
  neutral <- .fitted_random_effect_stan_data(
    draws_list = draw_list,
    stan_data = list(
      D = 1L,
      R_id = 1L,
      R_mk = 0L,
      Q_idm = 2L,
      indep_idmarker_cov = 0L
    ),
    reuse_fitted_re = FALSE
  )
  expect_equal(dim(neutral$z_L), c(2L, 3L))
})
