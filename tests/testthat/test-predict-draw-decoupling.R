test_that("prediction draw count resolves independently from posterior extraction", {
  expect_equal(.get_n_pred_draws(NULL, 12), 12)
  expect_equal(.get_n_pred_draws(7, 12), 7)
  expect_error(.get_n_pred_draws(0, 12), "must be >= 1")
  expect_error(.get_n_pred_draws(c(2, 3), 12), "single finite numeric")
})


test_that("prediction scale normalization supports multi-scale selection", {
  expect_equal(.normalize_prediction_scales("epred"), "epred")
  expect_equal(
    .normalize_prediction_scales(c("predict", "epred", "predict")),
    c("predict", "epred")
  )
  expect_error(.normalize_prediction_scales("bad_scale"), "should be one of")
})


test_that("prediction draw variable extraction keeps only requested scales", {
  standata_subject <- list(n_random_marker = 1L, n_random_marker_id = 0L)

  vars_epred <- .prediction_draw_variables("epred", standata_subject)
  expect_true("y_pred_epred" %in% vars_epred)
  expect_false("y_pred_linpred" %in% vars_epred)
  expect_false("y_pred" %in% vars_epred)
  expect_true("y_fit_epred" %in% vars_epred)
  expect_false("y_fit_linpred" %in% vars_epred)

  vars_predict <- .prediction_draw_variables("predict", standata_subject)
  expect_true("y_pred" %in% vars_predict)
  expect_false("y_fit_epred" %in% vars_predict)
  expect_false("y_fit_linpred" %in% vars_predict)

  vars_multi <- .prediction_draw_variables(c("linpred", "predict"), list(n_random_marker = 0L, n_random_marker_id = 1L))
  expect_true(all(c("y_pred_linpred", "y_pred", "y_fit_linpred", "z_w_lat", "z_L") %in% vars_multi))
  expect_false("y_fit_epred" %in% vars_multi)
})


test_that("prediction draw index supports downsample and upsample", {
  idx_down <- .prediction_draw_index(n_available = 10, n_target = 4, seed = 99)
  expect_length(idx_down, 4)
  expect_true(all(idx_down >= 1 & idx_down <= 10))
  expect_equal(length(unique(idx_down)), 4)

  idx_up <- .prediction_draw_index(n_available = 3, n_target = 8, seed = 99)
  expect_length(idx_up, 8)
  expect_true(all(idx_up >= 1 & idx_up <= 3))
})


test_that("draw-dependent containers are re-indexed to prediction draw count", {
  draws <- list(
    beta_fixed = matrix(seq_len(20), nrow = 10, ncol = 2),
    tau_id = matrix(seq_len(30), nrow = 10, ncol = 3),
    Lcorr_id = array(seq_len(40), dim = c(10, 2, 2)),
    lambda_vcov_reg = matrix(seq_len(20), nrow = 10, ncol = 2),
    bs_gamma_c = array(seq_len(30), dim = c(10, 1, 3))
  )

  idx <- c(2, 5, 9, 1)
  out <- .subset_draws_for_prediction(draws, idx)

  expect_equal(nrow(out$beta_fixed), length(idx))
  expect_equal(nrow(out$tau_id), length(idx))
  expect_equal(dim(out$Lcorr_id)[1], length(idx))
  expect_equal(nrow(out$lambda_vcov_reg), length(idx))
  expect_equal(dim(out$bs_gamma_c)[1], length(idx))

  expect_equal(out$beta_fixed[, 1], draws$beta_fixed[idx, 1])
  expect_equal(out$lambda_vcov_reg[, 1], draws$lambda_vcov_reg[idx, 1])
})


test_that("marker-id draw reconstruction retains the original-time covariance basis", {
  draws_matrix <- matrix(
    c(log(2), 3, log(4), 0, 0, 0, 5, 7),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "alpha_L[1]", "alpha_L[2]", "alpha_L[3]",
        "lambda_L[1]", "lambda_L[2]", "lambda_L[3]",
        "z_w_lat[1,1,1]", "z_w_lat[1,1,2]"
      )
    )
  )

  standata_subject <- list(
    n_random_marker_id = 2L,
    n_marker_types = 1L,
    alpha_vcov_reg = matrix(c(log(2), 3, log(4)), nrow = 1),
    beta_vcov_sd_flat = matrix(numeric(0), nrow = 1, ncol = 0),
    beta_vcov_corr_flat = matrix(numeric(0), nrow = 1, ncol = 0),
    lambda_vcov_reg = matrix(0, nrow = 1, ncol = 3),
    vec_cov_vcov_sd = numeric(0),
    vec_cov_vcov_corr = numeric(0),
    tau_marker = matrix(numeric(0), nrow = 1, ncol = 0),
    Lcorr_marker = array(numeric(0), dim = c(1, 0, 0)),
    B_cross = array(numeric(0), dim = c(1, 2, 0)),
    idx_row_cov = c(1L, 2L, 2L),
    idx_col_cov = c(1L, 1L, 2L),
    n_random_marker = 0L,
    flag_indep_marker_re = 1L,
    flag_allow_marker_crosscorr = 0L,
    vcov_diag_link = 1L,
    zidm_cols = c("(Intercept)", "time")
  )

  out <- .reconstruct_subject_marker_id_draws(draws_matrix, standata_subject, n_draws_target = 1)

  k21 <- tanh(3)
  l_i <- matrix(
    c(
      2, 0,
      4 * k21, 4 * sqrt(1 - k21^2)
    ),
    nrow = 2,
    byrow = TRUE
  )
  l_i_eff <- l_i
  expected_w <- as.numeric(l_i_eff %*% c(5, 7))

  expect_equal(unname(out$matrix[1, ]), expected_w, tolerance = 1e-8)
  expect_equal(out$corr[1, , ], l_i %*% t(l_i), tolerance = 1e-8)
})
