test_that("prediction draw count resolves independently from posterior extraction", {
  expect_equal(.resolve_n_pred_draws(NULL, 12), 12)
  expect_equal(.resolve_n_pred_draws(7, 12), 7)
  expect_error(.resolve_n_pred_draws(0, 12), "must be >= 1")
  expect_error(.resolve_n_pred_draws(c(2, 3), 12), "single finite numeric")
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
    tau_vcov_reg = seq_len(10),
    bs_gamma_c = array(seq_len(30), dim = c(10, 1, 3))
  )

  idx <- c(2, 5, 9, 1)
  out <- .subset_draws_for_prediction(draws, idx)

  expect_equal(nrow(out$beta_fixed), length(idx))
  expect_equal(nrow(out$tau_id), length(idx))
  expect_equal(dim(out$Lcorr_id)[1], length(idx))
  expect_equal(length(out$tau_vcov_reg), length(idx))
  expect_equal(dim(out$bs_gamma_c)[1], length(idx))

  expect_equal(out$beta_fixed[, 1], draws$beta_fixed[idx, 1])
  expect_equal(out$tau_vcov_reg, draws$tau_vcov_reg[idx])
})
