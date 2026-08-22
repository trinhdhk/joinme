test_that("supplied population parameters are fixed exactly and retained in truth", {
  simulation <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    n_id = 4,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.2),
    time_cens = 0.3,
    assoc = c("cv_mean", "cv_marker"),
    truth = jm_truth(
      longitudinal = list(
        intercept = c("(Intercept)" = 1.5),
        slope = c(time = -0.2, x1 = 0.4)
      ),
      survival = list(slope = c(x1 = 0.3, x2 = -0.1)),
      assoc_coef = list(slope = c(cv_mean = -0.7, cv_marker = 1)),
      marker_weights = list(intercept = 1.25, family = "normal"),
      re_params = list(
        id = list(sd = c(0.8, 0.4), corr = matrix(c(1, 0.2, 0.2, 1), 2)),
        marker = list(sd = 0.6, corr = matrix(1, 1, 1))
      )
    ),
    seed = 7331,
    use_mirai = FALSE
  )

  expect_equal(
    simulation$truth$beta_long[c("(Intercept)", "time", "x1")],
    c("(Intercept)" = 1.5, time = -0.2, x1 = 0.4)
  )
  expect_equal(simulation$truth$beta_event[c("x1", "x2")], c(x1 = 0.3, x2 = -0.1))
  expect_equal(
    simulation$truth$assoc_coefs[c("cv_mean", "cv_marker")],
    c(cv_mean = -0.7, cv_marker = 1)
  )
  expect_identical(simulation$truth$marker_weight_mean, c(shared = 1.25))
  expect_equal(simulation$truth$re_structure$id$sd, c(0.8, 0.4))
  expect_equal(simulation$truth$re_structure$marker$sd, 0.6)
})

test_that("omitted population parameters are drawn once and stored in truth", {
  simulation_arguments <- list(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.2),
    time_cens = 0.3,
    assoc = "cv_marker",
    truth = jm_truth(
      marker_weights = list(family = "normal"),
      basehaz = list(type = "weibull")
    ),
    use_mirai = FALSE
  )
  simulation_a <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7332)))
  simulation_b <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7333)))

  expect_true(all(is.finite(simulation_a$truth$beta_long)))
  expect_true(all(is.finite(simulation_a$truth$beta_event)))
  expect_true(is.finite(simulation_a$truth$assoc_coefs[["cv_marker"]]))
  expect_true(is.finite(simulation_a$truth$marker_weight_mean[["shared"]]))
  expect_true(all(is.finite(simulation_a$truth$re_structure$id$sd)))
  expect_true(all(is.finite(simulation_a$truth$baseline_hazard$parameters$shape)))
  expect_false(identical(
    simulation_a$truth$marker_weight_mean,
    simulation_b$truth$marker_weight_mean
  ))
  expect_false(identical(
    simulation_a$truth$beta_long,
    simulation_b$truth$beta_long
  ))
})

test_that("numeric marker-weight locations belong to truth rather than fitting priors", {
  expect_error(
    jm_priors(marker_weights = list(intercept = 1, family = "normal")),
    "fixed numerical value"
  )

  simulation <- simulate_joinme(
    n_id = 2,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.1),
    time_cens = 0.2,
    assoc = "cv_marker",
    truth = jm_truth(assoc_coef = list(slope = c(cv_marker = 0.4)),
      marker_weights = list(intercept = 1, family = "normal")
    ),
    seed = 7335,
    use_mirai = FALSE
  )
  expect_equal(simulation$truth$assoc_coefs[["cv_marker"]], 0.4)
  expect_equal(simulation$truth$marker_weight_mean[["shared"]], 1)
  expect_equal(
    simulation$truth$recovery$arguments$priors$marker_weights$intercept$mu,
    0
  )
  expect_equal(
    simulation$truth$recovery$arguments$priors$marker_weights$intercept$scale,
    2
  )
})

test_that("coefficient prior families generate one reproducible population vector", {
  population_distributions <- list(
    normal = prior_normal(mu = 0.5, scale = 0.2),
    student_t = prior_student_t(df = 7, mu = 0.5, scale = 0.2),
    laplace = prior_laplace(mu = 0.5, scale = 0.2),
    horseshoe = prior_horseshoe(global_scale = 0.25)
  )

  for (distribution_name in names(population_distributions)) {
    simulation_arguments <- list(
      formulaLong = y ~ 1 + time + (1 + time | id),
      formulaEvent = NULL,
      n_id = 3,
      families = "gaussian",
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      truth = jm_truth(longitudinal = list(
        intercept = population_distributions[[distribution_name]],
        slope = population_distributions[[distribution_name]]
      )),
      use_mirai = FALSE
    )
    first <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7340)))
    repeated <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7340)))
    independent <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7341)))

    expect_equal(first$truth$beta_long, repeated$truth$beta_long, info = distribution_name)
    expect_false(
      isTRUE(all.equal(first$truth$beta_long, independent$truth$beta_long)),
      info = distribution_name
    )
    expect_true(all(is.finite(first$truth$beta_long)), info = distribution_name)
  }
})

test_that("a bare component prior aligns coefficient-specific hyperparameters in model-matrix order", {
  population_prior <- prior_normal(
    mu = c(1, -0.2, 0.4),
    scale = c(0.1, 0.2, 0.3)
  )
  simulation <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = NULL,
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.1),
    truth = jm_truth(longitudinal = population_prior),
    seed = 7342,
    use_mirai = FALSE
  )
  prepared <- do.call(
    joinme,
    c(simulation$truth$recovery$arguments, list(fit = FALSE))
  )

  expect_equal(
    prepared$stan_data$prior_regression_mu[seq_along(simulation$truth$beta_long)],
    c(1, -0.2, 0.4)
  )
  expect_equal(
    prepared$stan_data$prior_regression_scale[seq_along(simulation$truth$beta_long)],
    c(0.1, 0.2, 0.3)
  )

  horseshoe_simulation <- simulate_joinme(
    formulaLong = simulation$truth$formulaLong,
    formulaEvent = NULL,
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.1),
    truth = jm_truth(longitudinal = prior_horseshoe(global_scale = 0.2)),
    seed = 7343,
    use_mirai = FALSE
  )
  horseshoe_prepared <- do.call(
    joinme,
    c(horseshoe_simulation$truth$recovery$arguments, list(fit = FALSE))
  )
  beta_horseshoe_groups <- horseshoe_prepared$stan_data$prior_regression_horseshoe_group_index[
    seq_along(horseshoe_simulation$truth$beta_long)
  ]
  expect_true(all(beta_horseshoe_groups > 0L))
  expect_length(unique(beta_horseshoe_groups), 1L)
})
