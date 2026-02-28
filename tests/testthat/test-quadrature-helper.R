test_that("gk_quadrature defaults to exact GK-15", {
  q <- gk_quadrature()
  expect_equal(q$n_gk, 15L)
  expect_equal(length(q$nodes), 15L)
  expect_equal(length(q$weights), 15L)
  expect_equal(sum(q$weights), 1, tolerance = 1e-12)
  expect_equal(q$rule, "gk15")
  expect_equal(q$panels, 1L)
})

test_that("gk_quadrature supports exact Boost-style rules", {
  q7 <- gk_quadrature(nodes = 7)
  expect_equal(q7$n_gk, 7L)
  expect_equal(length(q7$nodes), 7L)
  expect_equal(length(q7$weights), 7L)
  expect_equal(q7$rule, "gk7")
  expect_equal(sum(q7$weights), 1, tolerance = 1e-12)

  q <- gk_quadrature(nodes = 31)
  expect_equal(q$n_gk, 31L)
  expect_equal(q$panels, 1L)
  expect_equal(length(q$nodes), 31L)
  expect_equal(length(q$weights), 31L)
  expect_equal(q$rule, "gk31")
  expect_equal(sum(q$weights), 1, tolerance = 1e-12)
})

test_that("gk_quadrature rejects unsupported node requests", {
  expect_error(gk_quadrature(nodes = 32), "Allowed values")
  expect_error(gk_quadrature(nodes = 60), "Allowed values")
})
