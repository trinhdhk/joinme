test_that(".median_survival_time_by_draw extrapolates when threshold not reached on-grid", {
  s_mat <- rbind(
    c(0.95, 0.90, 0.85),
    c(0.80, 0.70, 0.60)
  )
  t_grid <- c(1, 2, 3)

  out <- joinme:::.median_survival_time_by_draw(
    survival_draw_matrix = s_mat,
    time_grid = t_grid,
    threshold = 0.5
  )

  expect_equal(length(out), 2)
  expect_true(all(is.finite(out)))
  expect_true(all(out > max(t_grid)))
})
