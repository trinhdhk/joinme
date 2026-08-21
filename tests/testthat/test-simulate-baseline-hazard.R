test_that("simulate_joinme supports piecewise baseline hazard", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    seed = 9001,
    formulaEvent = survival::Surv(time, event) ~ 1,
    assoc = c("cv_total"),
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      basehaz = list(type = "piecewise", breaks = c(2, 4), rates = c(0.05, 0.20, 0.40))
    )
  )

  bh <- sim$helpers$baseline_hazard
  vals <- bh(c(1, 3, 5))
  expect_equal(as.numeric(vals), c(0.05, 0.20, 0.40), tolerance = 1e-8)
  expect_true(all(vals > 0))
  expect_equal(sim$truth$baseline_hazard$type, "piecewise")
  expect_equal(sim$truth$baseline_hazard$parameters$rates, c(0.05, 0.20, 0.40))
})

test_that("simulate_joinme exposes the fitted-scale truth for a constant baseline", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    times_obs = c(0, 1),
    time_cens = 1,
    assoc = "cv_total",
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      basehaz = list(type = "constant", lambda = 0.1)
    ),
    seed = 9000,
    use_mirai = FALSE
  )

  expect_equal(sim$truth$baseline_hazard$parameters$rate, 0.1)
  expect_equal(
    unname(sim$truth$stan_fit$bs_gamma_c),
    log(0.1),
    tolerance = 1e-12
  )
})

test_that("simulate_joinme supports formula-based nonlinear baseline hazard", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    seed = 9002,
    formulaEvent = survival::Surv(time, event) ~ 1,
    assoc = c("cv_total"),
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      baseline = c("(Intercept)" = -2.0, "time" = -0.2, "I(time^2)" = 0.25),
      basehaz = ~ 1 + time + I(time^2)
    )
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

test_that("simulate_joinme spline baseline uses fixed coefficients from truth", {
  run_case <- function(kind) {
    sim <- simulate_joinme(
      n_id = 4,
      families = rep("student_t", 2),
      times_obs = seq(0, 8, length.out = 10),
      seed = 9010,
      formulaEvent = survival::Surv(time, event) ~ 1,
      assoc = c("cv_total"),
      truth = jm_truth(
        assoc_coef = list(slope = c(cv_total = 0)),
        baseline = -2,
        basehaz = list(type = kind, knots = c(3))
      ),
      integration_control = list(rel.tol = 1e-5, subdivisions = 10000L, stop.on.error = TRUE)
    )

    bh <- sim$helpers$baseline_hazard
    grid <- seq(0.2, 6, length.out = 12)
    v1 <- as.numeric(bh(grid))
    v2 <- as.numeric(bh(grid))

    expect_true(all(is.finite(v1)))
    expect_true(all(v1 > 0))
    expect_equal(v1, v2, tolerance = 1e-12)

    ch <- vapply(c(0.5, 1, 2, 4), function(tt) sim$helpers$cumhaz(1, tt), numeric(1))
    expect_true(all(is.finite(ch)))
    expect_true(all(diff(ch) > 0))
  }

  run_case("ns")
  run_case("bs")
})

test_that("simulate_joinme uses a formula baseline declared in jm_truth", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    times_obs = seq(0, 6, length.out = 8),
    seed = 9011,
    formulaEvent = survival::Surv(time, event) ~ 1,
    assoc = c("cv_total"),
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      baseline = c("(Intercept)" = -2.0, "time" = 0.3),
      basehaz = ~ 1 + time
    )
  )

  bh <- sim$helpers$baseline_hazard
  grid <- c(0.5, 1.5, 3.0)
  got <- as.numeric(bh(grid))
  expected <- exp(-2.0 + 0.3 * grid)
  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("simulate_joinme piecewise baselines use their rate parameters", {
  sim_a <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    seed = 9012,
    formulaEvent = survival::Surv(time, event) ~ 1,
    assoc = c("cv_total"),
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      basehaz = list(type = "piecewise", breaks = c(2, 4), rates = c(0.05, 0.20, 0.40))
    )
  )

  sim_b <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    seed = 9012,
    formulaEvent = survival::Surv(time, event) ~ 1,
    assoc = c("cv_total"),
    truth = jm_truth(
      assoc_coef = list(slope = c(cv_total = 0)),
      basehaz = list(type = "piecewise", breaks = c(2, 4), rates = c(0.05, 0.20, 0.40))
    )
  )

  grid <- c(1, 3, 5)
  expect_equal(
    as.numeric(sim_a$helpers$baseline_hazard(grid)),
    as.numeric(sim_b$helpers$baseline_hazard(grid)),
    tolerance = 1e-12
  )
})

test_that("simulate_joinme integrate path is stable across baseline hazard modes", {
  run_case <- function(mode, baseline_spec = NULL, formula_bh = NULL, beta_bh = NULL) {
    sim <- simulate_joinme(
      n_id = 6,
      families = rep("gaussian", 2),
      times_obs = seq(0, 6, length.out = 8),
      seed = 9013,
      formulaEvent = survival::Surv(time, event) ~ 1,
      assoc = c("cv_total"),
      truth = jm_truth(
        assoc_coef = list(slope = c(cv_total = 0)),
        baseline = beta_bh,
        basehaz = if (!is.null(formula_bh)) formula_bh else baseline_spec
      ),
      integration_control = list(rel.tol = 1e-5, subdivisions = 10000L, stop.on.error = TRUE)
    )

    bh <- sim$helpers$baseline_hazard
    grid <- c(0.2, 1, 2.5, 4.5)
    v1 <- as.numeric(bh(grid))
    v2 <- as.numeric(bh(grid))
    expect_true(all(is.finite(v1)))
    expect_true(all(v1 > 0))
    expect_equal(v1, v2, tolerance = 1e-12)

    ch <- vapply(c(0.5, 1, 2, 4), function(tt) sim$helpers$cumhaz(1, tt), numeric(1))
    expect_true(all(is.finite(ch)))
    expect_true(all(diff(ch) > 0))
  }

  run_case("constant", baseline_spec = list(type = "constant", rate = 0.08))
  run_case("linear", baseline_spec = list(type = "linear", intercept = -2.3, slope = 0.35))
  run_case("piecewise", baseline_spec = list(type = "piecewise", breaks = c(2, 4), rates = c(0.05, 0.2, 0.4)))
  run_case("weibull", baseline_spec = list(type = "weibull", shape = 1.4, scale = 6))
  run_case("formula", formula_bh = ~ 1 + time + I(time^2), beta_bh = c("(Intercept)" = -2.0, "time" = -0.1, "I(time^2)" = 0.05))
})
