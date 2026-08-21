test_that(".extract_matrix_from_stan pairs draw indices with posterior rows", {
  cn <- c(
    "y_pred_epred[1,1]", "y_pred_epred[2,1]",
    "y_pred_epred[1,2]", "y_pred_epred[2,2]"
  )

  draws_mat <- matrix(
    c(
      1, 3, 5, 7,
      2, 4, 6, 8
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, cn)
  )

  out <- .extract_matrix_from_stan(
    draws_mat = draws_mat,
    var_name = "y_pred_epred",
    N_cols = 2,
    N_rows = 2
  )

  expect_equal(dim(out), c(2, 2))
  expect_equal(out[, 1], c(1, 4))
  expect_equal(out[, 2], c(5, 8))
})

test_that(".extract_matrix_from_stan supports swapped index order", {
  cn <- c(
    "surv_prob[1,1]", "surv_prob[1,2]", "surv_prob[1,3]",
    "surv_prob[2,1]", "surv_prob[2,2]", "surv_prob[2,3]"
  )

  draws_mat <- matrix(
    c(
      10, 20, 30, 40, 50, 60,
      14, 24, 34, 44, 54, 64
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, cn)
  )

  out <- .extract_matrix_from_stan(
    draws_mat = draws_mat,
    var_name = "surv_prob",
    N_cols = 2,
    N_rows = 3
  )

  expect_equal(dim(out), c(3, 2))
  expect_equal(out[, 1], c(10, 24, 30))
  expect_equal(out[, 2], c(40, 54, 60))
})
