test_that("jm_prior accepts distributions and rejects fixed coefficients", {
  fitting_prior <- jm_prior(
    intercept = prior_normal(0, 2),
    slope = prior_student_t(5, 0, 1),
    assoc = list(slope = prior_laplace(0, 1))
  )

  expect_s3_class(fitting_prior, "joinme_priors")
  expect_error(jm_prior(longitudinal = c(0, 1)), "fixed numerical value")
  expect_error(
    jm_prior(marker_weights = list(intercept = 1, family = "normal")),
    "fixed numerical value"
  )
  expect_error(
    jm_prior(sigma = list(intercept = log(0.5))),
    "fixed numerical value"
  )
  expect_error(
    jm_truth(marker = 1),
    "centred unit-scale random-effect family"
  )
})

test_that("jm_truth retains fixed values and generating distributions", {
  declaration <- jm_truth(
    longitudinal = list(intercept = 1, slope = prior_normal(0.5, 0.1)),
    survival = c(treatment = -0.4),
    assoc_coef = c(cv_mean = 0.3),
    basehaz = list(type = "constant", rate = 0.08),
    re_params = list(id = list(sd = NULL, corr = NULL)),
    lkj = prior_lkj(2)
  )

  expect_s3_class(declaration, "joinme_truth")
  expect_identical(declaration$longitudinal$intercept, 1)
  expect_s3_class(declaration$longitudinal$slope, "joinme_prior_spec")
  expect_identical(declaration$survival, c(treatment = -0.4))
  expect_identical(declaration$assoc_coef, c(cv_mean = 0.3))
  expect_identical(declaration$basehaz$type, "constant")
  expect_equal(declaration$lkj$eta, 2)
})

test_that("truth distributions are drawn once per simulated data set", {
  generating_truth <- jm_truth(
    longitudinal = list(
      intercept = prior_normal(1, 0.2),
      slope = prior_normal(-0.3, 0.1)
    ),
    re_params = list(id = list(sd = NULL, corr = NULL)),
    lkj = prior_lkj(3)
  )
  arguments <- list(
    formulaLong = y ~ 1 + time + (1 + time | id),
    formulaEvent = NULL,
    n_id = 4,
    families = "gaussian",
    times_obs = c(0, 1),
    truth = generating_truth,
    use_mirai = FALSE
  )

  first <- do.call(simulate_joinme, c(arguments, list(seed = 8831)))
  repeated <- do.call(simulate_joinme, c(arguments, list(seed = 8831)))
  independent <- do.call(simulate_joinme, c(arguments, list(seed = 8832)))

  expect_equal(first$truth$beta_long, repeated$truth$beta_long)
  expect_equal(first$truth$re_structure$id$sd, repeated$truth$re_structure$id$sd)
  expect_equal(first$truth$re_structure$id$corr, repeated$truth$re_structure$id$corr)
  expect_false(isTRUE(all.equal(first$truth$beta_long, independent$truth$beta_long)))
  expect_true(all(first$truth$re_structure$id$sd > 0))
  expect_equal(diag(first$truth$re_structure$id$corr), c(1, 1))
})

test_that("simulation truth contains the complete baseline and random-effect declaration", {
  simulation <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id),
    formulaEvent = NULL,
    n_id = 3,
    families = "gaussian",
    times_obs = c(0, 1),
    truth = jm_truth(
      longitudinal = c(0.5, -0.2),
      basehaz = list(type = "weibull", shape = 1.3, scale = 7),
      re_params = list(id = list(sd = c(0.4, 0.2), corr = diag(2)))
    ),
    use_mirai = FALSE,
    seed = 8833
  )

  expect_equal(simulation$truth$beta_long, c("(Intercept)" = 0.5, time = -0.2))
  expect_equal(simulation$truth$re_structure$id$sd, c(0.4, 0.2))
  expect_identical(simulation$truth$declaration$basehaz$type, "weibull")
  expect_identical(simulation$truth$declaration$re_params$id$sd, c(0.4, 0.2))
})

test_that("simulation entry points expose one truth argument", {
  ordinary_arguments <- names(formals(simulate_joinme))
  mixture_arguments <- names(formals(simulate_joinme_mix))
  removed_arguments <- c("priors", "re_params", "h0", "baseline_hazard", "formulaBasehaz")

  expect_true("truth" %in% ordinary_arguments)
  expect_true("truth" %in% mixture_arguments)
  expect_length(intersect(ordinary_arguments, removed_arguments), 0L)
  expect_length(intersect(mixture_arguments, removed_arguments), 0L)
})
