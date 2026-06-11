test_that("simulate_joinme supports piecewise baseline hazard", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    seed = 9001,
    formulaEvent = survival::Surv(time, event) ~ 1,
    beta_event = c("(Intercept)" = 0),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0),
    baseline_hazard = list(type = "piecewise", breaks = c(2, 4), rates = c(0.05, 0.20, 0.40))
  )

  bh <- sim$helpers$baseline_hazard
  vals <- bh(c(1, 3, 5))
  expect_equal(as.numeric(vals), c(0.05, 0.20, 0.40), tolerance = 1e-8)
  expect_true(all(vals > 0))
})

test_that("simulate_joinme supports formula-based nonlinear baseline hazard", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    seed = 9002,
    formulaEvent = survival::Surv(time, event) ~ 1,
    beta_event = c("(Intercept)" = 0),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0),
    formulaBasehaz = ~ 1 + time + I(time^2),
    beta_basehaz = c("(Intercept)" = -2.0, "time" = -0.2, "I(time^2)" = 0.25)
  )

  bh <- sim$helpers$baseline_hazard
  vals <- as.numeric(bh(c(0.5, 1.0, 2.0, 3.0)))
  expect_true(all(vals > 0))
  expect_gt(vals[4], vals[2])

  ch <- vapply(c(0.5, 1.0, 1.5, 2.0), function(tt) sim$helpers$cumhaz(1, tt), numeric(1))
  dch <- diff(ch)
  expect_true(all(dch > 0))
  expect_gt(dch[length(dch)], dch[1])
})
