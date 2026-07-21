test_that("resolve_time_grid uses default horizon from tmax", {
  grid <- joinme:::.resolve_time_grid(
    times = NULL,
    id = "1",
    t_cond = 2,
    tmax_val = 10,
    time_horizon = 10,
    default_n = 10,
    min_points = 1,
    kind = "longitudinal"
  )

  expect_equal(min(grid), 2)
  expect_equal(max(grid), 10)
  expect_equal(length(grid), 10)
})


test_that("resolve_time_grid honors explicit shorter time_horizon", {
  grid <- joinme:::.resolve_time_grid(
    times = NULL,
    id = "1",
    t_cond = 2,
    tmax_val = 10,
    time_horizon = 2,
    default_n = 8,
    min_points = 1,
    kind = "survival"
  )

  expect_equal(min(grid), 2)
  expect_equal(max(grid), 4)
  expect_equal(length(grid), 8)
})


test_that("resolve_time_grid truncates horizon to training support", {
  expect_warning(
    grid <- joinme:::.resolve_time_grid(
      times = NULL,
      id = "1",
      t_cond = 2,
      tmax_val = 10,
      time_horizon = 20,
      default_n = 10,
      min_points = 1,
      kind = "longitudinal"
    ),
    "exceeds training time support"
  )

  expect_equal(min(grid), 2)
  expect_equal(max(grid), 10)
})


test_that("resolve_time_grid filters explicit times by condition and horizon", {
  expect_warning(
    grid <- joinme:::.resolve_time_grid(
      times = c(1, 2, 4, 6, 8),
      id = "1",
      t_cond = 2,
      tmax_val = 10,
      time_horizon = 4,
      default_n = 10,
      min_points = 1,
      kind = "survival"
    ),
    "exceeds configured horizon/support"
  )

  expect_equal(grid, c(2, 4, 6))
})
