test_that("fitted pwlin uses knots and ordered simplex metadata", {
  increasing <- build_standata_transforms(list(
    cv_total = list(
      type = "pwlin",
      knots = c(-2, -0.5, 1, 3),
      direction = "increasing"
    )
  ))
  decreasing <- build_standata_transforms(list(
    cv_total = list(
      type = "pwlin",
      cutpoints = c(-2, -0.5, 1, 3),
      direction = "decreasing"
    )
  ))

  expect_identical(increasing$tf_mode_cv_tot, 3)
  expect_identical(decreasing$tf_mode_cv_tot, 7)
  expect_equal(increasing$knots_cv, c(-2, -0.5, 1, 3))
  expect_equal(decreasing$knots_cv, c(-2, -0.5, 1, 3))
  expect_identical(increasing$estimate_spline_cv, 1L)
  expect_identical(decreasing$estimate_spline_cv, 1L)
  expect_identical(increasing$n_free_spline_cv, 3L)
  expect_identical(decreasing$n_free_spline_cv, 3L)
  expect_equal(increasing$coeff_cv, seq(0, 1, length.out = 4))
  expect_equal(decreasing$coeff_cv, seq(0, -1, length.out = 4))
})

test_that("legacy fitted pwlin y values no longer fix the association", {
  legacy <- build_standata_transforms(list(
    corr = list(
      type = "pwlin",
      x = c(-2, -1, 0, 1, 2),
      y = c(0.2, 0.5, 1, 0.5, 0.2)
    )
  ), n_corr_components = 2L)

  expect_identical(legacy$estimate_spline_corr, 1L)
  expect_identical(legacy$n_free_spline_corr, 4L)
  expect_equal(legacy$coeff_corr[1, ], seq(0, 1, length.out = 5))
  expect_equal(legacy$coeff_corr[2, ], seq(0, 1, length.out = 5))
  expect_false(isTRUE(all.equal(legacy$coeff_corr[1, ], c(0.2, 0.5, 1, 0.5, 0.2))))
})

test_that("legacy fitted pwlin infers decreasing direction from endpoints", {
  legacy <- build_standata_transforms(list(
    vcov = list(
      type = "pwlin",
      x = c(-1, 0, 1),
      y = c(4, 2, 1)
    )
  ), n_vcov_components = 1L)

  expect_identical(legacy$tf_mode_vcov, 7)
  expect_equal(as.numeric(legacy$coeff_vcov[1, ]), c(0, -0.5, -1))
})

test_that("fitted pwlin validates knots, direction, and legacy lengths", {
  expect_error(
    build_standata_transforms(list(cv_total = list(type = "pwlin", direction = "increasing"))),
    "require.*knots"
  )
  expect_error(
    build_standata_transforms(list(cv_total = list(type = "pwlin", knots = c(0, 0, 1)))),
    "strictly increasing"
  )
  expect_error(
    build_standata_transforms(list(cv_total = list(type = "pwlin", knots = c(0, 1), direction = "sideways"))),
    "increasing.*decreasing"
  )
  expect_error(
    build_standata_transforms(list(cv_total = list(type = "pwlin", x = c(0, 1, 2), y = c(0, 1)))),
    "one value per knot"
  )
})

test_that("piecewise-linear interpolation basis is constant-tailed and efficient", {
  basis <- joinme:::.pwlin_interpolation_basis(
    x_grid = c(-2, -1, -0.5, 0, 1, 2),
    knots = c(-1, 0, 1)
  )

  expect_equal(rowSums(basis), rep(1, 6))
  expect_equal(as.numeric(basis %*% c(0, 0.25, 1)), c(0, 0, 0.125, 0.25, 1, 1))
  expect_true(all(rowSums(basis != 0) <= 2L))
})

test_that("piecewise-linear plotting uses posterior ordinates", {
  object <- list(
    stan_data = list(
      knots_cv = c(-1, 0, 1),
      n_coeff_cv = 3L,
      coeff_cv = c(0, 0.5, 1)
    ),
    config = list(transforms_spec = list(
      cv_total = list(type = "pwlin", knots = c(-1, 0, 1))
    ))
  )
  posterior_ordinates <- rbind(c(0, 0.2, 1), c(0, 0.8, 1))
  testthat::local_mocked_bindings(
    .transform_coeff_draws = function(...) posterior_ordinates,
    .package = "joinme"
  )

  out <- joinme:::.pwlin_transform_matrix(
    object,
    term_key = "cv_total",
    x_grid = c(-2, -0.5, 0.5, 2),
    n_draws = 2L
  )

  expect_equal(out, rbind(c(0, 0.1, 0.6, 1), c(0, 0.4, 0.9, 1)))
})

test_that("fit and dynamic Stan paths recognise ordered pwlin mode", {
  stan_helper <- function(...) testthat::test_path("..", "..", "inst", "stan", "helper", ...)
  fit_data <- paste(readLines(stan_helper("data", "fit_data.stan"), warn = FALSE), collapse = "\n")
  dyn_data <- paste(readLines(stan_helper("data", "dynpred_data.stan"), warn = FALSE), collapse = "\n")
  transform_code <- paste(readLines(stan_helper("functions", "composite_transform.stanfunctions"), warn = FALSE), collapse = "\n")
  parameter_code <- paste(readLines(stan_helper("parameters", "joinme_fit_common.stan"), warn = FALSE), collapse = "\n")
  transformed_code <- paste(readLines(stan_helper("transformed_parameters", "fit_scaling_and_effects.stan"), warn = FALSE), collapse = "\n")
  prediction_code <- paste(readLines(testthat::test_path("..", "..", "R", "predict.R"), warn = FALSE), collapse = "\n")

  expect_match(fit_data, "upper=7")
  expect_match(dyn_data, "upper=7")
  expect_match(transform_code, "mode == 3 \\|\\| mode == 7")
  expect_match(parameter_code, "simplex.*pwlin_simplex_cv")
  expect_match(transformed_code, "tf_mode_cv_tot == 7")
  expect_match(transformed_code, "delta = pwlin_simplex_cv")
  expect_match(transformed_code, "delta = softmax\\(z_spline_cv\\)")
  expect_match(prediction_code, "get_transform_coeff_draws\\(\"coeff_cv_eff\"")
})

test_that("Stan estimates ordered pwlin ordinates in a fitted joint model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")

  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 3, length.out = 6),
    assoc = "cv_total",
    assoc_coefs = c(cv_total = 0.3),
    seed = 1729
  )
  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = "cv_total",
    families = rep("gaussian", 2),
    transforms = joinme_tf(cv_total = list(
      type = "pwlin",
      knots = c(0, 0.25, 0.75, 1.5),
      direction = "increasing"
    )),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0,
      seed = 1729
    )
  )

  ordinate_names <- paste0("coeff_cv_eff[", 1:4, "]")
  ordinate_draws <- joinme:::.get_draws_matrix(fit$fit, variables = ordinate_names)
  expect_s3_class(fit, "JoiNMeFit")
  expect_true(all(apply(ordinate_draws, 1L, function(x) all(diff(x) >= 0))))
  expect_equal(as.numeric(ordinate_draws[, 1L]), rep(0, nrow(ordinate_draws)))
  expect_equal(as.numeric(ordinate_draws[, 4L]), rep(1, nrow(ordinate_draws)), tolerance = 1e-10)
  expect_true(is.data.frame(summary(fit)$tables$piecewise_ordinates))
})

test_that("time-dependent AUC uses the shared dynamic-risk path", {
  object <- structure(list(), class = "JoiNMeFit")
  testthat::local_mocked_bindings(
    .resolve_train_data = function(object, newdataLong, newdataEvent, purpose) {
      list(newdataLong = data.frame(id = 1:4), newdataEvent = data.frame(id = 1:4))
    },
    .time_varying_concordance_single = function(...) {
      data.frame(
        risk = c(0.8, 0.6, 0.4, 0.6),
        event_window = c(1L, 1L, 0L, 0L),
        event_time = c(1.5, 1.8, 3, 4),
        time_horizon = rep(2, 4)
      )
    },
    .package = "joinme"
  )

  out <- auc(object, time_start = 1, time_horizon = 2)
  expect_s3_class(out, "tvAUC_JoiNMeFit")
  expect_equal(out$auc, 0.875)
  expect_equal(out$n_cases, 2L)
  expect_equal(out$n_controls, 2L)
  expect_equal(out$n_pairs, 4L)
})

test_that("piecewise summary reports posterior relative log-hazard ordinates", {
  object <- list(stan_data = list(
    knots_cv = c(-1, 0, 1),
    tf_mode_cv_tot = 3
  ))
  plot_data <- list(
    term_map = data.frame(term = "cv_total", variable = "alpha_cv_total"),
    coeff_draws = list(cv_total = matrix(c(2, 4), ncol = 1))
  )
  testthat::local_mocked_bindings(
    .get_association_plot_data = function(object, seed) plot_data,
    .assoc_coeff_draws = function(object, term, data, seed) c(2, 4),
    .association_transform_matrix = function(...) {
      rbind(c(0, 0.25, 1), c(0, 0.75, 1))
    },
    .package = "joinme"
  )

  out <- joinme:::.summarise_pwlin_ordinates(
    object,
    transform_specs = list(cv_total = list(
      type = "pwlin",
      knots = c(-1, 0, 1),
      direction = "increasing"
    ))
  )

  expect_equal(out$Estimate, c(0, 1.75, 3))
  expect_equal(out$Hazard.Ratio, round(exp(c(0, 1.75, 3)), 3))
  expect_equal(out$direction, rep("increasing", 3))
})
