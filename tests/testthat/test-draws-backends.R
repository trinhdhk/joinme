test_that("RStan draws are converted to indexed draws-matrix columns", {
  raw <- array(
    seq_len(16),
    dim = c(2, 2, 4),
    dimnames = list(
      iterations = NULL,
      chains = c("chain:1", "chain:2"),
      parameters = c("a", "z_u[1]", "z_u[2]", "z_L[1,1]")
    )
  )
  extract_call <- NULL

  local_mocked_bindings(
    extract = function(...) {
      extract_call <<- list(...)
      raw
    },
    .package = "rstan"
  )

  fit <- structure(list(), class = "stanfit")
  out <- .get_draws_matrix(fit, variables = c("z_u", "z_L"))

  expect_equal(dim(out), c(4, 4))
  expect_equal(
    colnames(out),
    c("a", "z_u[1]", "z_u[2]", "z_L[1,1]")
  )
  expect_identical(extract_call$pars, c("z_u", "z_L"))
  expect_false(extract_call$permuted)
  expect_false(extract_call$inc_warmup)
})


test_that("RStan draw subsampling uses the unified posterior representation", {
  raw <- array(
    seq_len(12),
    dim = c(3, 2, 2),
    dimnames = list(NULL, NULL, c("theta[1]", "theta[2]"))
  )

  local_mocked_bindings(
    extract = function(...) raw,
    .package = "rstan"
  )

  fit <- structure(list(), class = "stanfit")
  out <- suppressMessages(.get_draws_matrix(fit, draws = 4, seed = 2026))

  expect_equal(dim(out), c(4, 2))
  expect_equal(colnames(out), c("theta[1]", "theta[2]"))
})
