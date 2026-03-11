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
