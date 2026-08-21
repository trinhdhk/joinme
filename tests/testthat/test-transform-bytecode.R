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

  bc_pow <- parse_transform_expr(~ pow(x, 2))
  expect_equal(bc_pow, bc_power)
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

test_that("earlier unary bytecode is normalised with implicit PUSH_X", {
  earlier_softplus <- .normalize_bytecode_program(bytecode = 24L, const_data = numeric(0))

  expect_equal(earlier_softplus$bytecode, c(0L, 24L))
  expect_equal(eval_bytecode_scalar(2.5, earlier_softplus$bytecode, earlier_softplus$const_data), softplus(2.5))
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

test_that("the portable bytecode registry is complete and stable", {
  expect_identical(
    unname(.bytecode_opcodes()),
    0:27
  )
  expect_identical(
    names(.bytecode_opcodes())[c(1L, 2L, 27L, 28L)],
    c("PUSH_X", "PUSH_CONST", "PHI", "INV_PHI")
  )
})

test_that("the R interpreter can be loaded without JoiNMe internals", {
  module_file <- testthat::test_path(
    "..",
    "..",
    "R",
    "bytecode-interpreter.R"
  )
  isolated_module <- new.env(parent = baseenv())
  sys.source(module_file, envir = isolated_module)

  programme <- c(
    isolated_module$.bytecode_opcodes()[["PUSH_X"]],
    isolated_module$.bytecode_opcodes()[["PUSH_CONST"]],
    isolated_module$.bytecode_opcodes()[["ADD"]],
    isolated_module$.bytecode_opcodes()[["SOFTPLUS"]]
  )
  observed <- isolated_module$eval_bytecode_scalar(
    x = 2,
    bytecode = programme,
    const_data = 3
  )
  expect_equal(observed, log1p(exp(5)))
})

test_that("Stan entry points use the neutral interpreter module", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  entry_points <- file.path(
    stan_root,
    c(
      "joinme_fit_threading.stan",
      "joinme_mix_fit_threading.stan",
      "joinme_dynpred_threading.stan",
      "joinme_mix_dynpred_threading.stan"
    )
  )
  sources <- vapply(
    entry_points,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )

  expect_true(all(grepl(
    "include/etc/bytecode/interpreter.stanfunctions",
    sources,
    fixed = TRUE
  )))
})
