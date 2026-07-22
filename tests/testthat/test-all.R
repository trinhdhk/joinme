
testthat::test_that("all", {
  skip_on_cran()
  set.seed(42)
  sim <- simulate_joinme(
    n_id = 400,
    families = rep("student_t", 3),
    times_obs = seq(0, 8, length.out = 20),
    seed = 42,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.5),
    fixed_marker_weights = FALSE,
    shrinkage = 2L,
    baseline_hazard = list(type = "constant", lambda = 0.1),
    # beta_basehaz = c(-1),
    integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~  x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = "student_t",
    basehaz = joinme_basehaz(type = "formula", formula = ~ 1),
    transforms = joinme_tf(cv_total = "identity"),
    fixed_marker_weights = FALSE,
    shrinkage = 2L,
    priors = joinme_priors(
      beta = list(scale = 3),
      alpha = list(scale = 3),
      iota = list(scale = 3)
    ),
    control = list(
      engine = "rstan",
      parallel_chains = 1,
      threads_per_chain = 2,
      iter_warmup = 1000,
      iter_sampling = 3000,
      refresh = 200,
      adapt_delta = 0.7,
      max_treedepth = 11,
      seed = 421
    )
  )

  expect_s3_class(fit, "JoiNMeFit")
  expect_equal(sim$truth$shrinkage, fit$stan_data$shrinkage)
  expect_false(sim$truth$fixed_marker_weights)
  expect_equal(unname(sim$truth$marker_weights_base), fit$stan_data$marker_weights)
  expect_equal(
    unname(sim$truth$marker_weights),
    unname(sim$truth$marker_weights_base + sim$truth$marker_weights_latent)
  )
  expect_equal(sim$truth$baseline_hazard$parameters$rate, 0.1)
  expect_equal(unname(sim$truth$stan_fit$bs_gamma_c), log(0.1), tolerance = 1e-12)
  sum_obj <- summary(fit)
  expect_true(!is.null(sum_obj$tables))

  pred <- posterior_predict(
    fit,
    newdataLong = sim$dataLong[id==1],
    newdataEvent = sim$dataEvent[id==1],
    time_start = max(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 1, length.out = 20),
    control = list(
      n_samples = 20,
      chains = 1,
      iter_warmup = 500,
      iter_sampling = 500,
      refresh = 0
    )
  )

  expect_s3_class(pred, "JoiNMeDynPred")
  p <- plot(pred, type = c("longitudinal", "survival"), combined = FALSE)
  expect_true(inherits(p, "ggplot") || is.list(p))
})
