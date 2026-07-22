test_that("parse_transform_expr supports new bytecode ops", {
  bc_sigmoid <- parse_transform_expr(~ sigmoid(x))
  expect_equal(bc_sigmoid$bytecode, c(0L, 9L))

  bc_expit <- parse_transform_expr(~ expit(x))
  expect_equal(bc_expit$bytecode, c(0L, 9L))

  bc_softmax <- parse_transform_expr(~ SoftMax(x))
  expect_equal(bc_softmax$bytecode, c(0L, 9L))

  bc_softmax_lower <- parse_transform_expr(~ softmax(x))
  expect_equal(bc_softmax_lower$bytecode, bc_softmax$bytecode)

  bc_softplus <- parse_transform_expr(~ softplus(x))
  expect_equal(bc_softplus$bytecode, c(0L, 24L))

  bc_log1p <- parse_transform_expr(~ log1p_exp(x))
  expect_equal(bc_log1p$bytecode, c(0L, 24L))

  bc_cbrt <- parse_transform_expr(~ cbrt(x))
  expect_equal(bc_cbrt$bytecode, c(0L, 25L))

  bc_power <- parse_transform_expr(~ power(x, 2))
  expect_equal(bc_power$bytecode, c(0L, 1L, 12L))
  expect_equal(bc_power$const_data, 2)
})

test_that("parse_transform_expr treats unary minus as 0-x", {
  bc_unary <- parse_transform_expr(~ -x)
  bc_binary <- parse_transform_expr(~ 0 - x)

  expect_equal(bc_unary$bytecode, c(1L, 0L, 3L))
  expect_equal(bc_unary$const_data, 0)
  expect_equal(bc_unary$bytecode, bc_binary$bytecode)
  expect_equal(bc_unary$const_data, bc_binary$const_data)

  expect_equal(eval_bytecode_scalar(2.5, bc_unary$bytecode, bc_unary$const_data), -2.5)
  expect_equal(
    eval_bytecode_vector(c(-2, 0, 3), bc_unary$bytecode, bc_unary$const_data),
    c(2, 0, -3)
  )
})

test_that("legacy unary bytecode is normalised with implicit PUSH_X", {
  legacy_softplus <- .normalize_bytecode_program(bytecode = 24L, const_data = numeric(0))

  expect_equal(legacy_softplus$bytecode, c(0L, 24L))
  expect_equal(eval_bytecode_scalar(2.5, legacy_softplus$bytecode, legacy_softplus$const_data), softplus(2.5))
  expect_equal(
    eval_bytecode_vector(c(-2, 0, 3), 24L, numeric(0)),
    softplus(c(-2, 0, 3))
  )
})

test_that("bytecode evaluator applies per-node affine shifts", {
  tf <- joinme_tf(
    cv_total = ~ softplus(expit(x, intercept = TRUE, slope = TRUE), intercept = TRUE, slope = TRUE)
  )
  bc <- parse_transform_expr(tf$cv_total$expr, iota_nodes = tf$cv_total$iota_nodes)

  val <- eval_bytecode_scalar(
    x = 0.3,
    bytecode = bc$bytecode,
    const_data = bc$const_data,
    iota_intercepts = c(0.4, -0.2),
    iota_slopes = c(1.5, 0.8),
    op_iota_intercept_idx = bc$op_iota_intercept_idx,
    op_iota_slope_idx = bc$op_iota_slope_idx
  )

  expected <- softplus(-0.2 + 0.8 * stats::plogis(0.4 + 1.5 * 0.3))
  expect_equal(val, expected)
})
