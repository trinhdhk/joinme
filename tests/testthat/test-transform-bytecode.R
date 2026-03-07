test_that("parse_transform_expr supports new bytecode ops", {
  bc_sigmoid <- parse_transform_expr(~ sigmoid(x))
  expect_equal(bc_sigmoid$opcodes, c(0L, 9L))
  expect_equal(bc_sigmoid$bytecode, bc_sigmoid$opcodes)

  bc_expit <- parse_transform_expr(~ expit(x))
  expect_equal(bc_expit$opcodes, c(0L, 9L))

  bc_softplus <- parse_transform_expr(~ softplus(x))
  expect_equal(bc_softplus$opcodes, c(0L, 24L))

  bc_log1p <- parse_transform_expr(~ log1p_exp(x))
  expect_equal(bc_log1p$opcodes, c(0L, 24L))

  bc_cbrt <- parse_transform_expr(~ cbrt(x))
  expect_equal(bc_cbrt$opcodes, c(0L, 25L))

  bc_power <- parse_transform_expr(~ power(x, 2))
  expect_equal(bc_power$opcodes, c(0L, 1L, 12L))
  expect_equal(bc_power$const_data, 2)
})

test_that("parse_transform_expr treats unary minus as 0-x", {
  bc_unary <- parse_transform_expr(~ -x)
  bc_binary <- parse_transform_expr(~ 0 - x)

  expect_equal(bc_unary$opcodes, c(1L, 0L, 3L))
  expect_equal(bc_unary$const_data, 0)
  expect_equal(bc_unary$opcodes, bc_binary$opcodes)
  expect_equal(bc_unary$const_data, bc_binary$const_data)

  expect_equal(eval_bytecode_scalar(2.5, bc_unary$bytecode, bc_unary$const_data), -2.5)
  expect_equal(
    eval_bytecode_vector(c(-2, 0, 3), bc_unary$bytecode, bc_unary$const_data),
    c(2, 0, -3)
  )
})
