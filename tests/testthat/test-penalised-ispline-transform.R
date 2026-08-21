test_that("penalised_ispline_transform uses anchored endpoint coefficients", {
  skip_if_not_installed("splines2")

  tf <- penalised_ispline_transform(
    x = c(0, 0.25, 0.5, 0.6, 0.8, 1),
    y = c(0, 0.5, 0.6, 0.75, 0.8, 1),
    knots = seq(0, 1, 0.2),
    degree = 2,
    lambda = 1.5
  )

  expect_identical(tf$type, "ispline")
  expect_equal(tf$knots, seq(0, 1, 0.2))
  expect_length(tf$coeff, length(tf$knots) + tf$degree - 1)
  expect_equal(tf$coeff[1], 0, tolerance = 1e-8)
  expect_equal(tf$coeff[length(tf$coeff)], 1, tolerance = 1e-8)
  expect_true(all(diff(tf$coeff) >= -1e-8))
})

test_that("expit-domain penalised spline keeps the same anchored monotone contract", {
  skip_if_not_installed("splines2")

  tf <- .make_penalised_ispline_transform(list(
    type = "ispline_expit_penalised",
    x = seq(0.05, 0.95, length.out = 8),
    y = seq(0.05, 0.95, length.out = 8)^0.8,
    knots = seq(0.05, 0.95, length.out = 5),
    degree = 2,
    lambda = 1
  ))

  expect_identical(tf$type, "ispline_expit")
  expect_true(all(tf$knots >= 0 & tf$knots <= 1))
  expect_equal(tf$coeff[1], 0, tolerance = 1e-8)
  expect_equal(tf$coeff[length(tf$coeff)], 1, tolerance = 1e-8)
  expect_true(all(diff(tf$coeff) >= -1e-8))
})

test_that("penalised spline supports decreasing monotone fits", {
  skip_if_not_installed("splines2")

  z_grid <- seq(-2, 2, length.out = 8)
  tf <- .make_penalised_ispline_transform(list(
    type = "ispline_expit_penalised",
    x = stats::plogis(z_grid),
    y = exp(-0.5 * z_grid),
    knots = seq(0.05, 0.95, length.out = 5),
    degree = 2,
    lambda = 1,
    direction = "decreasing"
  ))

  tf_fun <- .make_assoc_transform(tf, "vcov")
  vals <- tf_fun(z_grid)

  expect_identical(tf$direction, "decreasing")
  expect_identical(tf$spline_direction, -1L)
  expect_equal(tf$coeff[1], 0, tolerance = 1e-8)
  expect_equal(tf$coeff[length(tf$coeff)], 1, tolerance = 1e-8)
  expect_true(all(diff(tf$coeff) >= -1e-8))
  expect_true(all(diff(vals) <= 1e-8))
})

test_that("Stan-estimated penalised spline and expit variant share the same monotone setup", {
  tf_raw <- .make_stan_penalised_ispline_transform(list(
    type = "ispline_penalised",
    x = seq(-2, 2, length.out = 10),
    n_knots = 5,
    degree = 3,
    lambda = 0.5
  ))
  tf_expit <- .make_stan_penalised_ispline_transform(list(
    type = "ispline_expit_penalised",
    x = seq(0.05, 0.95, length.out = 10),
    n_knots = 5,
    degree = 3,
    lambda = 0.5
  ))

  expect_identical(tf_raw$type, "ispline")
  expect_identical(tf_expit$type, "ispline_expit")
  expect_equal(tf_raw$n_free_spline, length(tf_raw$coeff) - 1L)
  expect_equal(tf_expit$n_free_spline, length(tf_expit$coeff) - 1L)
  expect_true(all(tf_expit$knots >= 0 & tf_expit$knots <= 1))
})

test_that("Stan-estimated penalised spline stores decreasing direction metadata", {
  tf <- .make_stan_penalised_ispline_transform(list(
    type = "ispline_expit_penalised",
    x = seq(0.05, 0.95, length.out = 10),
    n_knots = 5,
    degree = 3,
    lambda = 0.5,
    direction = "decreasing"
  ))

  expect_identical(tf$type, "ispline_expit")
  expect_identical(tf$direction, "decreasing")
  expect_identical(tf$spline_direction, -1L)
})
