test_that("plot.JoiNMeFit uses fitted draw builder for fitted plots", {
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL

  testthat::local_mocked_bindings(
    longitudinal_plot = function(object, ...) {
      captured$args <- list(...)
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(1, 2, 1.5, 2.5)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  out <- plot(
    fit,
    type = "longitudinal",
    subject = c(1, 2),
    scale = "epred",
    longitudinal_times = c(0, 0.5, 1),
    longitudinal_points = 40,
    draws = 20,
    show_data = FALSE
  )

  expect_s3_class(out, "ggplot")
  expect_equal(captured$args$scale, "epred")
  expect_equal(captured$args$draws, 20)
  expect_equal(as.character(captured$args$subject), c("1", "2"))
  expect_equal(captured$args$longitudinal_times, c(0, 0.5, 1))
  expect_equal(captured$args$longitudinal_points, 40)
})

test_that("fitted sample builder aligns subset subjects to Stan longitudinal rows", {
  base_draws <- matrix(
    c(
      0.10, 0.80, 0.70,
      0.20, 0.75, 0.65,
      0.15, 0.78, 0.68,
      0.12, 0.82, 0.72,
      0.18, 0.79, 0.69
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("beta[1]", "surv_prob_event[1]", "surv_prob_event[2]"))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      posterior::as_draws_matrix(base_draws)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      if (is.null(variables)) {
        return(base_draws)
      }
      keep <- intersect(variables, colnames(base_draws))
      if (!length(keep)) {
        return(matrix(0, nrow = nrow(base_draws), ncol = 0))
      }
      base_draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      n_id = 2L,
      N = 4L,
      D = 1L,
      P = 1L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      X_obs = matrix(1, nrow = 4, ncol = 1),
      Z_id_obs = matrix(0, nrow = 4, ncol = 0),
      Z_mk_obs = matrix(0, nrow = 4, ncol = 0),
      Z_idm_obs = matrix(0, nrow = 4, ncol = 0),
      marker_levels = "m1",
      link_long = 1L,
      inv_link_n_ops = 0L,
      inv_link_n_const = 0L,
      inv_link_ops = matrix(0L, nrow = 1, ncol = 1),
      inv_link_const = matrix(0, nrow = 1, ncol = 1)
    ),
    formulaLong = y ~ 1,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(draws_default = 5),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1, formulaEvent = survival::Surv(time, event) ~ 1)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 2, 2),
      time = c(0, 1, 0, 1),
      marker = "m1",
      y = c(1.0, 1.1, 0.9, 1.0)
    ),
    dataEvent = data.frame(
      id = c(1, 2),
      time = c(2.0, 2.5),
      event = c(0L, 1L)
    )
  )

  pred <- .build_JoiNMefit_fitted_plot_samples(
    x = fit,
    which = "survival",
    subject = 1,
    scale = "epred",
    draws = 5,
    seed = 1,
    ci_levels = c(0.5, 0.95)
  )

  expect_s3_class(pred, "JoiNMeDynPred")
  expect_setequal(unique(as.character(pred$quantiles$survival$id)), "1")
  expect_null(pred$quantiles$longitudinal)
  expect_null(pred$quantiles$cumhaz)
  expect_true(all(as.character(unique(pred$data$longitudinal$id)) == "1"))
})

test_that("fitted sample builder returns non-degenerate survival time grid", {
  base_draws <- matrix(
    c(
      0.10, 0.80, 0.70,
      0.20, 0.75, 0.65,
      0.15, 0.78, 0.68,
      0.12, 0.82, 0.72,
      0.18, 0.79, 0.69
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("beta[1]", "surv_prob_event[1]", "surv_prob_event[2]"))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      posterior::as_draws_matrix(base_draws)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      if (is.null(variables)) {
        return(base_draws)
      }
      keep <- intersect(variables, colnames(base_draws))
      if (!length(keep)) {
        return(matrix(0, nrow = nrow(base_draws), ncol = 0))
      }
      base_draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      n_id = 2L,
      N = 4L,
      D = 1L,
      P = 1L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      X_obs = matrix(1, nrow = 4, ncol = 1),
      Z_id_obs = matrix(0, nrow = 4, ncol = 0),
      Z_mk_obs = matrix(0, nrow = 4, ncol = 0),
      Z_idm_obs = matrix(0, nrow = 4, ncol = 0),
      marker_levels = "m1",
      link_long = 1L,
      inv_link_n_ops = 0L,
      inv_link_n_const = 0L,
      inv_link_ops = matrix(0L, nrow = 1, ncol = 1),
      inv_link_const = matrix(0, nrow = 1, ncol = 1)
    ),
    formulaLong = y ~ 1,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(draws_default = 5),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1, formulaEvent = survival::Surv(time, event) ~ 1)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 2, 2),
      time = c(0, 1, 0, 1),
      marker = "m1",
      y = c(1.0, 1.1, 0.9, 1.0)
    ),
    dataEvent = data.frame(
      id = c(1, 2),
      time = c(2.0, 2.5),
      event = c(0L, 1L)
    )
  )

  pred <- .build_JoiNMefit_fitted_plot_samples(
    x = fit,
    which = "survival",
    subject = 1,
    scale = "epred",
    draws = 5,
    seed = 1,
    ci_levels = c(0.5, 0.95)
  )

  surv_id1 <- pred$quantiles$survival[pred$quantiles$survival$id == "1", , drop = FALSE]
  expect_gt(length(unique(surv_id1$time)), 2)
})

test_that("plot.JoiNMeFit can plot association curves", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      map <- data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      vals <- matrix(seq(0.2, 0.6, length.out = 8), ncol = 1)
      colnames(vals) <- "cv_total"
      if (!is.null(term)) {
        map <- map[map$term %in% term, , drop = FALSE]
      }
      list(draws = vals, term_map = map)
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 1, 1.5))
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "association", association_options = list(association_term = "cv_total"))
  expect_s3_class(p, "ggplot")
})

test_that("plot.JoiNMeFit supports all diagnostic plot types with parameter filters", {
  draws_array <- posterior::as_draws_array(array(
    data = c(
      seq(0.1, 1.2, length.out = 24),
      seq(1.1, 2.2, length.out = 24),
      seq(-0.5, 0.6, length.out = 24)
    ),
    dim = c(12, 2, 3),
    dimnames = list(NULL, NULL, c("alpha", "beta", "gamma"))
  ))

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      out <- draws_array
      if (!is.null(variables)) {
        out <- posterior::subset_draws(out, variable = intersect(variables, posterior::variables(out)))
      }
      if (!is.null(draws) && is.finite(draws) && !isTRUE(keep_chains)) {
        n_keep <- min(as.integer(draws), posterior::ndraws(out))
        out <- posterior::subset_draws(out, draw = seq_len(n_keep))
      }
      out
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  scalar_types <- c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd")
  for (diag_type in scalar_types) {
    p <- plot(fit, type = diag_type, pars = "alpha")
    expect_s3_class(p, "ggplot")
    expect_equal(unique(as.character(p$data$variable)), "alpha")
  }

  running_types <- c("running_mean", "running_quantile")
  for (diag_type in running_types) {
    p <- plot(
      fit,
      type = diag_type,
      regex_pars = "alpha|beta",
      max_vars = 2,
      quantile_probs = c(0.25, 0.5, 0.75),
      draws = 8
    )
    expect_s3_class(p, "ggplot")
    expect_true(all(unique(as.character(p$data$variable)) %in% c("alpha", "beta")))
  }
})

test_that("plot.JoiNMeFit filters fitted longitudinal plots by marker", {
  testthat::local_mocked_bindings(
    longitudinal_plot = function(object, ...) {
      args <- list(...)
      ggplot2::ggplot(
        data.frame(marker = as.character(args$marker %||% NA_character_), x = 1, y = 1),
        ggplot2::aes(x, y, color = marker)
      ) + ggplot2::geom_point()
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 1, 1),
      time = c(0, 1, 0, 1),
      marker = c("m1", "m1", "m2", "m2"),
      y = c(1, 2, 1.5, 2.5)
    ),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "longitudinal", subject = 1, marker = "m1", show_data = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(unique(as.character(p$data$marker)), "m1")

  p_all <- plot(fit, type = "longitudinal", subject = 1, marker = NA, show_data = FALSE)
  expect_true(all(is.na(p_all$data$marker)))
})

test_that("plot.JoiNMeFit routes single longitudinal requests through longitudinal_plot", {
  captured <- new.env(parent = emptyenv())
  captured$called <- FALSE

  testthat::local_mocked_bindings(
    longitudinal_plot = function(object, ...) {
      captured$called <- TRUE
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "longitudinal")
  expect_true(captured$called)
  expect_s3_class(p, "ggplot")
})

test_that("plot.JoiNMeFit forwards mcmc requests to mcmc_plot", {
  captured <- new.env(parent = emptyenv())
  captured$type <- NULL
  captured$variable <- NULL

  testthat::local_mocked_bindings(
    mcmc_plot = function(object, pars = NA, type = "intervals", variable = NULL, regex = FALSE,
                         fixed = FALSE, draws = NULL, seed = 1, ...) {
      captured$type <- type
      captured$variable <- variable
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "mcmc", pars = "time", mcmc_type = "trace")
  expect_identical(captured$type, "trace")
  expect_identical(captured$variable, "time")
  expect_s3_class(p, "ggplot")
})

test_that("plot.JoiNMeFit rejects removed conditioning and prediction arguments", {
  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + trt + time,
    formulaEvent = survival::Surv(time, event) ~ age,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + trt + time, formulaEvent = survival::Surv(time, event) ~ age)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1),
      time = c(0, 1),
      marker = factor(c("m1", "m1"), levels = c("m1", "m2")),
      trt = factor(c("A", "A"), levels = c("A", "B")),
      y = c(1, 2)
    ),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L, age = 55)
  )

  expect_error(plot(fit, type = "longitudinal", condition = list(trt = "B")), "unused argument")
  expect_error(plot(fit, type = "longitudinal", conditioning = "last"), "unused argument")
  expect_error(plot(fit, type = "longitudinal", times = seq(0, 1, length.out = 5)), "unused argument")
  expect_error(plot(fit, type = "longitudinal", time_horizon = 2), "unused argument")
  expect_error(plot(fit, type = "longitudinal", pred_control = list(n_samples = 10)), "unused argument")
})

test_that("fitted longitudinal trajectories use a smooth common time design and retain observed responses", {
  fitted_draws <- matrix(
    c(
      0.0, 0.8,
      0.1, 1.0,
      -0.1, 1.2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("beta[1]", "beta[2]"))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      posterior::as_draws_matrix(fitted_draws)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      retained_variables <- intersect(variables %||% colnames(fitted_draws), colnames(fitted_draws))
      fitted_draws[, retained_variables, drop = FALSE]
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      n_id = 2L,
      N = 4L,
      D = 1L,
      P = 2L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      X_obs = cbind(1, c(0, 1, 0.2, 0.8)),
      Z_id_obs = matrix(0, nrow = 4, ncol = 0),
      Z_mk_obs = matrix(0, nrow = 4, ncol = 0),
      Z_idm_obs = matrix(0, nrow = 4, ncol = 0),
      marker_levels = "m1",
      link_long = 1L,
      inv_link_n_ops = 0L,
      inv_link_n_const = 0L,
      inv_link_ops = matrix(0L, nrow = 1, ncol = 1),
      inv_link_const = matrix(0, nrow = 1, ncol = 1),
      tmax = 1
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(draws_default = 3),
    call = quote(joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 2, 2),
      time = c(0, 1, 0.2, 0.8),
      marker = "m1",
      y = c(10, 11, 20, 21)
    ),
    dataEvent = data.frame(id = c(1, 2), time = c(2, 2), event = 0L)
  )

  fitted_prediction <- .build_JoiNMefit_fitted_plot_samples(
    x = fit,
    which = "longitudinal",
    subject = NULL,
    scale = "epred",
    draws = 3,
    seed = 1,
    ci_levels = c(0.5, 0.95),
    longitudinal_times = c(0, 0.5, 1)
  )

  subject_times <- split(
    fitted_prediction$quantiles$longitudinal$time,
    fitted_prediction$quantiles$longitudinal$id
  )
  expect_true(all(vapply(subject_times, identical, logical(1), c(0, 0.5, 1))))
  expect_s3_class(fitted_prediction, "JoiNMeDynPred")
  expect_s3_class(fitted_prediction, "R6")
  expect_equal(fitted_prediction$metadata$marker_levels, "m1")
  expect_null(fitted_prediction$quantiles$survival)
  expect_null(fitted_prediction$quantiles$cumhaz)
  expect_equal(fitted_prediction$data$longitudinal$y, c(10, 11, 20, 21))

  subject_plot <- plot(
    fitted_prediction,
    type = "longitudinal",
    subject = 1,
    scale = "epred",
    smooth_trajectory = FALSE
  )
  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), subject_plot$layers)
  expect_length(point_layers, 1L)
  expect_equal(point_layers[[1L]]$data$value, c(10, 11))
})


test_that("longitudinal heatmap aggregates smooth trajectories without requiring identical marker availability", {
  common_times <- c(0, 0.5, 1)
  posterior_prediction <- list(
    draws = list(
      longitudinal = list(
        `1` = list(
          epred = list(
            matrix = rbind(c(0, 1, 2, 4, 5, 6), c(0, 2, 4, 4, 6, 8)),
            time = rep(common_times, 2),
            marker_idx = rep(1:2, each = 3),
            scale = "epred"
          )
        ),
        `2` = list(
          epred = list(
            matrix = rbind(c(2, 3, 4), c(2, 4, 6)),
            time = common_times,
            marker_idx = rep(1L, 3),
            scale = "epred"
          )
        )
      )
    ),
    data = list(
      longitudinal = data.frame(marker = factor(c("m1", "m2"), levels = c("m1", "m2")))
    ),
    metadata = list(marker_var = "marker")
  )

  heatmap_data <- .JoiNMefit_longitudinal_heatmap_data(
    posterior_prediction = posterior_prediction,
    prediction_scale = "epred",
    sign_threshold = 0.05,
    prediction_times = common_times
  )

  expect_setequal(as.character(unique(heatmap_data$marker)), c("m1", "m2"))
  expect_equal(sort(unique(heatmap_data$time)), common_times)
  expect_equal(heatmap_data$change[as.character(heatmap_data$marker) == "m1"], c(0, 1.5, 3))
  expect_true(all(heatmap_data$alpha > 0))

  testthat::local_mocked_bindings(
    .build_JoiNMefit_fitted_plot_samples = function(...) posterior_prediction,
    .package = "joinme"
  )
  heatmap_plot <- .plot_JoiNMefit_longitudinal_heatmap(
    fitted_model = NULL,
    subject_ids = NULL,
    marker_levels = NULL,
    prediction_scale = "epred",
    posterior_draws = 2,
    prediction_times = common_times,
    number_time_points = 3,
    sign_threshold = 0.05,
    random_seed = 1,
    credible_levels = 0.95,
    plot_theme = ggplot2::theme_bw
  )
  tile_layer <- ggplot2::ggplot_build(heatmap_plot)$data[[1L]]
  expect_equal(nrow(tile_layer), nrow(heatmap_data))
  expect_true(all(tile_layer$alpha >= 0.3))
})


test_that("JoiNMeDynPred plotting preserves custom data-variable metadata", {
  quantiles_longitudinal <- data.frame(
    id = "A",
    time = c(0, 1),
    marker = "m1",
    q2.5 = c(0, 0.5),
    q50 = c(0.2, 0.8),
    q97.5 = c(0.4, 1.1),
    mean = c(0.2, 0.8),
    scale = "epred"
  )
  prediction <- JoiNMeDynPred$new(
    predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
    quantiles = list(
      longitudinal = quantiles_longitudinal,
      longitudinal_fitted = NULL,
      longitudinal_marker_pop = NULL,
      longitudinal_overall_pop = NULL,
      survival = NULL,
      cumhaz = NULL
    ),
    draws = list(longitudinal = list(), longitudinal_fitted = list(), survival = list(), cumhaz = list()),
    data = list(
      longitudinal = data.frame(
        patient = c("A", "A", "B"),
        visit_time = c(0, 1, 0),
        outcome = c("m1", "m1", "m1"),
        response = c(10, 11, 99)
      ),
      event = data.frame(patient = c("A", "B"), event_time = 2, event = 0L)
    ),
    metadata = list(
      id_var = "patient",
      time_var = "visit_time",
      marker_var = "outcome",
      marker_levels = "m1",
      response_var = "response",
      scale = "epred",
      scales = "epred",
      ci_levels = 0.95,
      conditioning_time = 0,
      conditioning_time_by_id = c(A = 0),
      source = "fit_samples"
    ),
    call = quote(predict(fit)),
    tmax = 1,
    n_samples = 2
  )

  subject_plot <- plot(prediction, type = "longitudinal", subject = "A", smooth_trajectory = FALSE)
  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), subject_plot$layers)

  expect_length(point_layers, 1L)
  expect_equal(point_layers[[1L]]$data$value, c(10, 11))
})


test_that("event-only JoiNMeDynPred objects do not require a longitudinal scale", {
  quantiles_survival <- data.frame(
    id = 1,
    time = c(0, 1),
    q2.5 = c(1, 0.6),
    q50 = c(1, 0.8),
    q97.5 = c(1, 0.95),
    mean = c(1, 0.8)
  )
  prediction <- JoiNMeDynPred$new(
    predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
    quantiles = list(longitudinal = NULL, survival = quantiles_survival, cumhaz = NULL),
    draws = list(longitudinal = NULL, survival = list(), cumhaz = list()),
    data = list(longitudinal = NULL, event = data.frame(id = 1, time = 1, event = 0L)),
    metadata = list(ci_levels = 0.95, id_var = "id", time_var = "time", marker_var = "marker"),
    call = quote(predict(fit)),
    tmax = 1,
    n_samples = 2
  )

  survival_plot_object <- plot(prediction, type = "survival")
  expect_s3_class(survival_plot_object, "ggplot")
})

test_that("plot.JoiNMeFit longitudinal heatmap style orders markers and attenuates nonsignificant changes", {
  testthat::local_mocked_bindings(
    longitudinal_plot = function(object, ...) {
      args <- list(...)
      dat <- data.frame(time = c(0, 1), marker = c("m1", "m2"), change = c(0.5, 0.1), alpha = c(1, 0.3))
      p <- ggplot2::ggplot(dat, ggplot2::aes(time, marker, fill = change, alpha = alpha)) + ggplot2::geom_tile()
      attr(p, "longitudinal_style") <- args$longitudinal_style
      p
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 1, 2, 2, 2),
      time = c(0, 1, 0, 0, 1, 0),
      marker = factor(c("m1", "m1", "m2", "m1", "m1", "m2"), levels = c("m1", "m2")),
      y = c(1, 2, 1.2, 1.1, 2.1, 1.3)
    ),
    dataEvent = data.frame(id = c(1, 2), time = c(2.5, 2.8), event = c(0L, 1L))
  )

  p <- plot(
    fit,
    type = "longitudinal",
    longitudinal_style = "heatmap",
    scale = "epred",
    threshold = 0.05
  )

  expect_s3_class(p, "ggplot")
  expect_equal(attr(p, "longitudinal_style"), "heatmap")
  expect_true(any(p$data$alpha == 1))
  expect_true(any(p$data$alpha == 0.3))
  expect_true(all(p$data$alpha > 0))
  expect_true(all(c("m1", "m2") %in% unique(as.character(p$data$marker))))
})

test_that("plot.JoiNMeFit no longer accepts condition for heatmap plotting", {
  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + trt + time,
    formulaEvent = survival::Surv(time, event) ~ trt,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + trt + time, formulaEvent = survival::Surv(time, event) ~ trt)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1),
      time = c(0, 1),
      marker = factor(c("m1", "m2"), levels = c("m1", "m2")),
      trt = factor(c("A", "A"), levels = c("A", "B")),
      y = c(1, 1.2)
    ),
    dataEvent = data.frame(id = 1, time = 2, event = 0L, trt = factor("A", levels = c("A", "B")) )
  )

  cond <- data.frame(
    trt = factor(c("A", "B"), levels = c("A", "B")),
    row.names = c("Arm A", "Arm B")
  )

  expect_error(
    plot(fit, type = "longitudinal", longitudinal_style = "heatmap", scale = "epred", condition = cond, threshold = 0.05),
    "unused argument"
  )
})

test_that("plot.JoiNMeFit keeps longitudinal_heatmap as a compatibility alias", {
  testthat::local_mocked_bindings(
    longitudinal_plot = function(object, ...) {
      args <- list(...)
      dat <- data.frame(time = c(0, 1), marker = c("m1", "m1"), change = c(0.0, 0.5), alpha = c(1, 1))
      p <- ggplot2::ggplot(dat, ggplot2::aes(time, marker, fill = change, alpha = alpha)) + ggplot2::geom_tile()
      attr(p, "longitudinal_style") <- args$longitudinal_style
      p
    },
    .package = "joinme"
  )

  fit <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = 1, time = c(0, 1), marker = factor(c("m1", "m1"), levels = "m1"), y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 2, event = 0L)
  )

  p <- plot(fit, type = "longitudinal_heatmap", scale = "epred")
  expect_s3_class(p, "ggplot")
  expect_equal(attr(p, "longitudinal_style"), "heatmap")
})

test_that("make_conditions mirrors brms output", {
  testthat::skip_if_not_installed("brms")

  dat <- data.frame(
    trt = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
    age = c(50, 60, 55, 65)
  )

  expect_equal(
    JoiNMe::make_conditions(dat, vars = "trt"),
    brms::make_conditions(dat, vars = "trt")
  )
})

test_that("plot.JoiNMeFit association curves respect marker selection", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      map <- data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      vals <- matrix(seq(0.2, 0.6, length.out = 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(draws = vals, term_map = map)
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 1, 1.5))
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 1, 1),
      time = c(0, 1, 0, 1),
      marker = c("m1", "m1", "m2", "m2"),
      y = c(1, 2, 1.5, 2.5)
    ),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p_one <- plot(fit, type = "association", association_options = list(association_term = "cv_total"), marker = "m1")
  expect_s3_class(p_one, "ggplot")
  expect_equal(unique(as.character(p_one$data$marker)), "m1")

  expect_no_warning({
    p_all <- plot(fit, type = "association", association_options = list(association_term = "cv_total"), marker = NA)
  })
  expect_s3_class(p_all, "ggplot")
  expect_setequal(unique(as.character(p_all$data$marker)), c("m1", "m2"))
})

test_that("plot.JoiNMeFit rejects expit spline knots outside [0, 1]", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      map <- data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      vals <- matrix(seq(0.2, 0.6, length.out = 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(draws = vals, term_map = map)
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_expit", knots = c(-0.1, 0.5, 1.1), coeff = c(0, 0.4, 0.8), degree = 1)
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  expect_error(
    plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "transform")),
    "must lie on the expit scale"
  )
})

test_that("plot.JoiNMeFit weighted cv_total curves are marker-specific", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(
      marker_levels = c("m1", "m2"),
      marker_weights = c(1, 2)
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x - 0.5))
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 1, 1),
      time = c(0, 1, 0, 1),
      marker = c("m1", "m1", "m2", "m2"),
      y = c(1, 2, 1.5, 2.5)
    ),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- c(0.5)
  p <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "cv_total",
      association_metric = "transform",
      association_grid = x_grid
    ),
    marker = NA
  )

  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)
  vals <- stats::setNames(p$data[[median_col]], as.character(p$data$marker))
  expect_equal(unname(vals[["m2"]]), 2 * unname(vals[["m1"]]), tolerance = 1e-8)
})

test_that("plot.JoiNMeFit cv_total default grid uses observed y support", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x - 0.5))
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(10, 12, 20, 22)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  p <- plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "transform"))
  expect_gt(max(p$data$x), 5)
  expect_gt(min(p$data$x), 5)
})

test_that("plot.JoiNMeFit cv_total warns and falls back to knot support when y exceeds spline support", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_penalised", knots = c(0, 0.5, 1), degree = 2)
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(10, 12, 20, 22)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  expect_warning({
    p <- plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "transform"))
  }, "extends beyond the fitted spline knot range")
  expect_gte(min(p$data$x), 0)
  expect_lte(max(p$data$x), 1)
})

test_that("plot.JoiNMeFit association curves respect nonlinear functional transforms", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "corr"
      list(
        draws = vals,
        term_map = data.frame(term = "corr", variable = "alpha_corr", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        corr = list(expr = ~ expit(-x), type = "functional")
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-1, 1, length.out = 5)
  p <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "corr",
      association_metric = "transform",
      association_grid = x_grid
    )
  )

  expected <- 1 / (1 + exp(x_grid))
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)
  expect_equal(p$data[[median_col]], expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(diff(p$data[[median_col]]), rep(diff(p$data[[median_col]])[1], length(diff(p$data[[median_col]]))))))
})

test_that("plot.JoiNMeFit association_options range controls the raw evaluation grid", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x))
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "cv_total",
      association_metric = "transform",
      association_range = c(-4, 4),
      association_points = 5
    )
  )

  x_grid <- seq(-4, 4, length.out = 5)
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)

  expect_equal(p$data$x, x_grid, tolerance = 1e-8)
  expect_equal(p$data[[median_col]], stats::plogis(x_grid), tolerance = 1e-8)
})

test_that("plot.JoiNMeFit rejects simultaneous association grid and range", {
  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "identity")
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  expect_error(
    plot(
      fit,
      type = "association",
      association_options = list(
        association_term = "cv_total",
        association_grid = c(-1, 0, 1),
        association_range = c(-2, 2)
      )
    ),
    "association_grid.*association_range|association_range.*association_grid"
  )
})

test_that("plot.JoiNMeFit association curves use fitted spline coefficients", {
  skip_if_not_installed("splines2")

  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .get_draws_matrix = function(fit, variables, draws = NULL, seed = 1, ...) {
      out <- matrix(rep(c(0, 1, 1), each = 8), nrow = 8)
      colnames(out) <- variables
      out
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(
      n_coeff_cv = 3L,
      knots_cv = c(-1, 0, 1),
      spline_degree_cv = 1L
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_penalised")
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-1, 1, length.out = 5)
  p <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "cv_total",
      association_metric = "transform",
      association_grid = x_grid
    )
  )

  basis <- splines2::iSpline(
    x_grid,
    knots = 0,
    degree = 1,
    intercept = TRUE,
    Boundary.knots = c(-1, 1)
  )
  expected <- as.numeric(as.matrix(basis) %*% c(0, 1, 1))
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)
  expect_equal(p$data[[median_col]], expected, tolerance = 1e-8)
})

test_that("plot.JoiNMeFit association curves evaluate expit splines on the transformed basis scale", {
  skip_if_not_installed("splines2")

  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "corr"
      list(
        draws = vals,
        term_map = data.frame(term = "corr", variable = "alpha_corr", stringsAsFactors = FALSE)
      )
    },
    .get_draws_matrix = function(fit, variables, draws = NULL, seed = 1, ...) {
      out <- matrix(rep(c(0, 1, 1), each = 8), nrow = 8)
      colnames(out) <- variables
      out
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(
      n_coeff_corr = 3L,
      knots_corr = c(0.2, 0.5, 0.8),
      spline_degree_corr = 1L
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        corr = list(type = "ispline_expit_penalised")
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-2, 2, length.out = 5)
  p <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "corr",
      association_metric = "transform",
      association_grid = x_grid
    )
  )

  basis_x <- pmin(pmax(stats::plogis(x_grid), 0.2), 0.8)
  basis <- splines2::iSpline(
    basis_x,
    knots = 0.5,
    degree = 1,
    intercept = TRUE,
    Boundary.knots = c(0.2, 0.8)
  )
  expected <- as.numeric(as.matrix(basis) %*% c(0, 1, 1))
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)

  expect_equal(p$data$x, x_grid, tolerance = 1e-8)
  expect_equal(p$data[[median_col]], expected, tolerance = 1e-8)
})

test_that("plot.JoiNMeFit expands vcov association plots to all components", {
  payload <- list(
    coeff_draws = list(
      vcov = structure(
        matrix(c(rep(1, 6), rep(2, 6)), nrow = 6, ncol = 2),
        dimnames = list(NULL, c("alpha_vcov_eff[1]", "alpha_vcov_eff[2]"))
      )
    ),
    marker_weight_draws = NULL,
    transform_coeff_draws = list(),
    transform_specs = list(),
    support = data.frame(
      term = c("vcov[1]", "vcov[2]"),
      marker = c("all", "all"),
      lower = c(-1, 0.1),
      upper = c(1, 1.5),
      source = c("theoretical", "model_implied"),
      stringsAsFactors = FALSE
    ),
    term_map = data.frame(
      term = c("vcov[1]", "vcov[2]"),
      variable = c("alpha_vcov_eff[1]", "alpha_vcov_eff[2]"),
      stringsAsFactors = FALSE
    )
  )

  testthat::local_mocked_bindings(
    .get_association_plot_data = function(x, seed = 1) payload,
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      assoc = list(vcov = TRUE),
      transforms = list(),
      transforms_spec = list(vcov = list(type = "identity"))
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-1, 1, length.out = 4)
  plots <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "vcov",
      association_metric = "hazard",
      association_grid = x_grid
    ),
    combined = FALSE
  )

  expect_type(plots, "list")
  expect_setequal(names(plots), c("vcov[1]", "vcov[2]"))

  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)
  expect_equal(plots[["vcov[1]"]]$data[[median_col]], x_grid, tolerance = 1e-8)
  expect_equal(plots[["vcov[2]"]]$data[[median_col]], 2 * x_grid, tolerance = 1e-8)
})

test_that("plot.JoiNMeFit uses component-specific vcov support ranges", {
  payload <- list(
    coeff_draws = list(
      vcov = structure(
        matrix(c(rep(-1, 6), rep(-2, 6)), nrow = 6, ncol = 2),
        dimnames = list(NULL, c("alpha_vcov_eff[1]", "alpha_vcov_eff[2]"))
      )
    ),
    marker_weight_draws = NULL,
    transform_coeff_draws = list(),
    transform_specs = list(vcov = list(type = "identity")),
    support = data.frame(
      term = c("vcov[1]", "vcov[2]"),
      marker = c("all", "all"),
      lower = c(-0.8, 0.2),
      upper = c(0.6, 1.1),
      source = c("model_implied", "model_implied"),
      stringsAsFactors = FALSE
    ),
    term_map = data.frame(
      term = c("vcov[1]", "vcov[2]"),
      variable = c("alpha_vcov_eff[1]", "alpha_vcov_eff[2]"),
      stringsAsFactors = FALSE
    )
  )

  testthat::local_mocked_bindings(
    .get_association_plot_data = function(x, seed = 1) payload,
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      assoc = list(vcov = TRUE),
      transforms = list(),
      transforms_spec = list(vcov = list(type = "identity"))
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  plots <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "vcov",
      association_metric = "hazard"
    ),
    combined = FALSE
  )

  expect_equal(range(plots[["vcov[1]"]]$data$x), c(-0.8, 0.6), tolerance = 1e-8)
  expect_equal(range(plots[["vcov[2]"]]$data$x), c(0.2, 1.1), tolerance = 1e-8)
})

test_that("plot.JoiNMeFit zero-references covariance-style hazard contributions", {
  payload <- list(
    coeff_draws = list(
      vcov = structure(
        matrix(rep(1, 6), nrow = 6, ncol = 1),
        dimnames = list(NULL, "alpha_vcov_eff[1]")
      )
    ),
    marker_weight_draws = NULL,
    transform_coeff_draws = list(),
    transform_specs = list(vcov = list(type = "functional", expr = ~ exp(-x))),
    support = data.frame(
      term = "vcov",
      marker = "all",
      lower = 0,
      upper = 1,
      source = "theoretical",
      stringsAsFactors = FALSE
    ),
    term_map = data.frame(
      term = "vcov[1]",
      variable = "alpha_vcov_eff[1]",
      stringsAsFactors = FALSE
    )
  )

  testthat::local_mocked_bindings(
    .get_association_plot_data = function(x, seed = 1) payload,
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      assoc = list(vcov = TRUE),
      transforms = list(),
      transforms_spec = list(vcov = list(type = "functional", expr = ~ exp(-x)))
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(0, 1, length.out = 4)
  p_hazard <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "vcov",
      association_metric = "hazard",
      association_grid = x_grid
    ),
    combined = FALSE
  )
  if (is.list(p_hazard) && !inherits(p_hazard, "ggplot")) {
    p_hazard <- p_hazard[[1]]
  }

  p_transform <- plot(
    fit,
    type = "association",
    association_options = list(
      association_term = "vcov",
      association_metric = "transform",
      association_grid = x_grid
    ),
    combined = FALSE
  )
  if (is.list(p_transform) && !inherits(p_transform, "ggplot")) {
    p_transform <- p_transform[[1]]
  }

  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)
  expect_equal(p_transform$data[[median_col]], exp(-x_grid), tolerance = 1e-8)
  expect_equal(p_hazard$data[[median_col]], exp(-x_grid) - 1, tolerance = 1e-8)
})

test_that("association plot payload prefers effective covariance coefficients", {
  draws_obj <- posterior::as_draws_matrix(stats::setNames(
    data.frame(
      raw = c(0.1, 0.1, 0.1),
      eff = c(0.4, 0.5, 0.6),
      check.names = FALSE
    ),
    c("alpha_vcov[1]", "alpha_vcov_eff[1]")
  ))

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      as.matrix(posterior::subset_draws(draws_obj, variable = variables))
    },
    .model_implied_support = function(fit, stan_data, config, dataLong, seed = 1) {
      data.frame(term = "vcov", marker = "all", lower = -1, upper = 1, source = "mock", stringsAsFactors = FALSE)
    },
    .package = "joinme"
  )

  payload <- JoiNMe:::.build_association_plot_data(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(Q_idm = 1L, indep_idmarker_cov = 1L, D = 1L),
    config = list(
      assoc = list(vcov = TRUE),
      transforms_spec = list(vcov = list(type = "identity"))
    ),
    dataLong = data.frame(id = 1L, marker = "m1", time = 0, y = 0),
    seed = 1
  )

  expect_equal(payload$term_map$term, "vcov[1]")
  expect_equal(payload$term_map$variable, "alpha_vcov_eff[1]")
  expect_equal(as.numeric(payload$coeff_draws$vcov[, 1]), c(0.4, 0.5, 0.6))
})

test_that("plot.JoiNMeFit defaults expit-spline fallback grids on the raw scale", {
  testthat::local_mocked_bindings(
    extract.JoiNMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      vals <- matrix(rep(1, 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(
        draws = vals,
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    },
    .package = "joinme"
  )

  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_expit_penalised", knots = c(0.2, 0.5, 0.8), degree = 1)
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(10, 12, 20, 22)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  expect_warning({
    p <- plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "transform"))
  }, "extends beyond the fitted spline knot range")

  expect_lt(min(p$data$x), 0)
  expect_gt(max(p$data$x), 0)
  expect_equal(range(p$data$x), stats::qlogis(c(0.2, 0.8)), tolerance = 1e-6)
})

test_that("plot.JoiNMeFit uses cached association plotting payload when available", {
  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x))
      ),
      association_plot_data = list(
        coeff_draws = list(cv_total = matrix(rep(2, 6), ncol = 1)),
        marker_weight_draws = matrix(rep(1, 6), ncol = 1),
        transform_coeff_draws = list(),
        transform_specs = list(cv_total = list(type = "functional", expr = ~ expit(x))),
        support = data.frame(
          term = "cv_total",
          marker = "m1",
          lower = -3,
          upper = 3,
          source = "model_implied",
          stringsAsFactors = FALSE
        ),
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "hazard"))
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)

  expect_equal(range(p$data$x), c(-3, 3), tolerance = 1e-8)
  expect_equal(p$data[[median_col]], 2 * (1 / (1 + exp(-p$data$x))), tolerance = 1e-8)
})

test_that("plot.JoiNMeFit applies iota affine shifts to functional association transforms", {
  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x), estimate_iota_intercept = TRUE, estimate_iota_slope = TRUE)
      ),
      association_plot_data = list(
        coeff_draws = list(cv_total = matrix(rep(1, 6), ncol = 1)),
        marker_weight_draws = matrix(rep(1, 6), ncol = 1),
        transform_coeff_draws = list(),
        transform_iota_draws = list(cv_total = list(intercept = rep(0.5, 6), slope = rep(2, 6))),
        transform_specs = list(cv_total = list(type = "functional", expr = ~ expit(x), estimate_iota_intercept = TRUE, estimate_iota_slope = TRUE)),
        support = data.frame(
          term = "cv_total",
          marker = "m1",
          lower = -1,
          upper = 1,
          source = "model_implied",
          stringsAsFactors = FALSE
        ),
        term_map = data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      )
    ),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "association", association_options = list(association_term = "cv_total", association_metric = "transform"))
  median_col <- JoiNMe:::.quantile_name_from_prob(0.5)

  expect_equal(p$data[[median_col]], stats::plogis(0.5 + 2 * p$data$x), tolerance = 1e-8)
})
