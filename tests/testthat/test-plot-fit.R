test_that("plot.JoinMeFit reuses fitted-data prediction plotting", {
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL

  testthat::local_mocked_bindings(
    predict.JoinMeFit = function(object, newdataLong, newdataEvent, process, pred_type, scale,
                                 times, time_start, time_horizon, ci_levels, control, seed, ...) {
      captured$args <- list(
        process = process,
        scale = scale,
        time_start = time_start,
        ids = unique(newdataEvent$id)
      )

      time_grid <- seq(0, 1, length.out = 5)
      quant_long <- do.call(rbind, lapply(unique(newdataEvent$id), function(id) {
        data.frame(
          id = id,
          time = time_grid,
          marker = "m1",
          scale = "epred",
          q2.5 = seq(0.1, 0.3, length.out = length(time_grid)),
          q50 = seq(0.2, 0.4, length.out = length(time_grid)),
          q97.5 = seq(0.3, 0.5, length.out = length(time_grid)),
          mean = seq(0.2, 0.4, length.out = length(time_grid)),
          sd = rep(0.05, length(time_grid)),
          stringsAsFactors = FALSE
        )
      }))
      quant_surv <- do.call(rbind, lapply(unique(newdataEvent$id), function(id) {
        data.frame(
          id = id,
          time = time_grid,
          q2.5 = pmax(0, 0.8 - time_grid),
          q50 = pmax(0, 0.9 - time_grid),
          q97.5 = pmax(0, 1.0 - time_grid),
          mean = pmax(0, 0.9 - time_grid),
          sd = rep(0.03, length(time_grid)),
          stringsAsFactors = FALSE
        )
      }))

      joinme::JoinMeDynPred$new(
        predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
        quantiles = list(
          longitudinal = quant_long,
          longitudinal_fitted = quant_long,
          longitudinal_marker_pop = NULL,
          longitudinal_overall_pop = NULL,
          survival = quant_surv,
          cumhaz = NULL
        ),
        draws = list(longitudinal = list(), longitudinal_fitted = list(), survival = list(), cumhaz = list()),
        data = list(longitudinal = newdataLong, event = newdataEvent),
        metadata = list(
          scales = "epred",
          scale = "epred",
          ci_levels = ci_levels,
          conditioning_time = 0.5,
          conditioning_time_by_id = stats::setNames(rep(0.5, length(unique(newdataEvent$id))), unique(newdataEvent$id)),
          id_var = "id",
          time_var = "time",
          marker_var = "marker",
          response_var = "y"
        ),
        call = quote(predict(fit_obj)),
        tmax = 1,
        n_samples = if (is.null(control$n_samples)) 10 else control$n_samples
      )
    },
    .package = "joinme"
  )

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(transforms = list()),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(1, 2, 1.5, 2.5)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  out <- plot(
    fit,
    type = c("longitudinal", "survival"),
    subject = c(1, 2),
    combined = FALSE,
    show_data = FALSE,
    pred_control = list(n_samples = 20)
  )

  expect_true(is.list(out) || inherits(out, "gg"))
  expect_equal(sort(captured$args$process), c("event", "longitudinal"))
  expect_equal(captured$args$scale, c("epred", "linpred", "predict"))
  expect_equal(as.numeric(captured$args$time_start[c("1", "2")]), c(1, 1))
})

test_that("plot.JoinMeFit can plot association curves", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 1, 1.5))
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "association", association_term = "cv_total")
  expect_s3_class(p, "ggplot")
})

test_that("plot.JoinMeFit filters fitted longitudinal plots by marker", {
  captured <- new.env(parent = emptyenv())
  captured$process <- NULL

  testthat::local_mocked_bindings(
    predict.JoinMeFit = function(object, newdataLong, newdataEvent, process, pred_type, scale,
                                 times, time_start, time_horizon, ci_levels, control, seed, ...) {
      captured$process <- process
      time_grid <- seq(0, 1, length.out = 5)
      quant_long <- do.call(rbind, lapply(c("m1", "m2"), function(marker) {
        data.frame(
          id = 1,
          time = time_grid,
          marker = marker,
          scale = "epred",
          q2.5 = seq(0.1, 0.3, length.out = length(time_grid)),
          q50 = seq(0.2, 0.4, length.out = length(time_grid)),
          q97.5 = seq(0.3, 0.5, length.out = length(time_grid)),
          mean = seq(0.2, 0.4, length.out = length(time_grid)),
          sd = rep(0.05, length(time_grid)),
          stringsAsFactors = FALSE
        )
      }))

      joinme::JoinMeDynPred$new(
        predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
        quantiles = list(
          longitudinal = quant_long,
          longitudinal_fitted = quant_long,
          longitudinal_marker_pop = NULL,
          longitudinal_overall_pop = NULL,
          survival = NULL,
          cumhaz = NULL
        ),
        draws = list(longitudinal = list(), longitudinal_fitted = list(), survival = list(), cumhaz = list()),
        data = list(longitudinal = newdataLong, event = newdataEvent),
        metadata = list(
          scales = "epred",
          scale = "epred",
          ci_levels = ci_levels,
          conditioning_time = 0.5,
          conditioning_time_by_id = c(`1` = 0.5),
          id_var = "id",
          time_var = "time",
          marker_var = "marker",
          response_var = "y"
        ),
        call = quote(predict(fit_obj)),
        tmax = 1,
        n_samples = 10
      )
    },
    .package = "joinme"
  )

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(transforms = list(), transforms_spec = list()),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
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
  expect_equal(sort(captured$process), c("event", "longitudinal"))

  p_all <- plot(fit, type = "longitudinal", subject = 1, marker = NA, show_data = FALSE)
  expect_setequal(unique(as.character(p_all$data$marker)), c("m1", "m2"))
})

test_that("plot.JoinMeFit association curves respect marker selection", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                                 term = NULL, variable = NULL, draws = NULL, seed = 1, keep_chains = TRUE, ...) {
      map <- data.frame(term = "cv_total", variable = "alpha_cv_total", stringsAsFactors = FALSE)
      vals <- matrix(seq(0.2, 0.6, length.out = 8), ncol = 1)
      colnames(vals) <- "cv_total"
      list(draws = vals, term_map = map)
    },
    .package = "joinme"
  )

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 1, 1.5))
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(
      id = c(1, 1, 1, 1),
      time = c(0, 1, 0, 1),
      marker = c("m1", "m1", "m2", "m2"),
      y = c(1, 2, 1.5, 2.5)
    ),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p_one <- plot(fit, type = "association", association_term = "cv_total", marker = "m1")
  expect_s3_class(p_one, "ggplot")
  expect_equal(unique(as.character(p_one$data$marker)), "m1")

  expect_no_warning({
    p_all <- plot(fit, type = "association", association_term = "cv_total", marker = NA)
  })
  expect_s3_class(p_all, "ggplot")
  expect_setequal(unique(as.character(p_all$data$marker)), c("m1", "m2"))
})

test_that("plot.JoinMeFit weighted cv_total curves are marker-specific", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(
      marker_levels = c("m1", "m2"),
      marker_weights = c(1, 2)
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x - 0.5))
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
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
    association_term = "cv_total",
    association_metric = "transform",
    association_grid = x_grid,
    marker = NA
  )

  median_col <- joinme:::.quantile_name_from_prob(0.5)
  vals <- stats::setNames(p$data[[median_col]], as.character(p$data$marker))
  expect_equal(unname(vals[["m2"]]), 2 * unname(vals[["m1"]]), tolerance = 1e-8)
})

test_that("plot.JoinMeFit cv_total default grid uses observed y support", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x - 0.5))
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(10, 12, 20, 22)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  p <- plot(fit, type = "association", association_term = "cv_total", association_metric = "transform")
  expect_gt(max(p$data$x), 5)
  expect_gt(min(p$data$x), 5)
})

test_that("plot.JoinMeFit cv_total warns and falls back to knot support when y exceeds spline support", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_penalised", knots = c(0, 0.5, 1), degree = 2)
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1, 2, 2), time = c(0, 1, 0, 1), marker = "m1", y = c(10, 12, 20, 22)),
    dataEvent = data.frame(id = c(1, 2), time = c(1.2, 1.3), event = c(0L, 1L))
  )

  expect_warning({
    p <- plot(fit, type = "association", association_term = "cv_total", association_metric = "transform")
  }, "extends beyond the fitted spline knot range")
  expect_gte(min(p$data$x), 0)
  expect_lte(max(p$data$x), 1)
})

test_that("plot.JoinMeFit association curves respect nonlinear functional transforms", {
  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        corr = list(expr = ~ expit(-x), type = "functional")
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-1, 1, length.out = 5)
  p <- plot(
    fit,
    type = "association",
    association_term = "corr",
    association_metric = "transform",
    association_grid = x_grid
  )

  expected <- 1 / (1 + exp(x_grid))
  median_col <- joinme:::.quantile_name_from_prob(0.5)
  expect_equal(p$data[[median_col]], expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(diff(p$data[[median_col]]), rep(diff(p$data[[median_col]])[1], length(diff(p$data[[median_col]]))))))
})

test_that("plot.JoinMeFit association curves use fitted spline coefficients", {
  skip_if_not_installed("splines2")

  testthat::local_mocked_bindings(
    extract.JoinMeFit = function(object, what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
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

  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(
      n_coeff_cv = 3L,
      knots_cv = c(-1, 0, 1),
      spline_degree_cv = 1L
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "ispline_penalised")
      )
    ),
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  x_grid <- seq(-1, 1, length.out = 5)
  p <- plot(
    fit,
    type = "association",
    association_term = "cv_total",
    association_metric = "transform",
    association_grid = x_grid
  )

  basis <- splines2::iSpline(
    x_grid,
    knots = 0,
    degree = 1,
    intercept = TRUE,
    Boundary.knots = c(-1, 1)
  )
  expected <- as.numeric(as.matrix(basis) %*% c(0, 1, 1))
  median_col <- joinme:::.quantile_name_from_prob(0.5)
  expect_equal(p$data[[median_col]], expected, tolerance = 1e-8)
})

test_that("plot.JoinMeFit uses cached association plotting payload when available", {
  fit <- joinme::JoinMeFit$new(
    fit = NULL,
    stan_data = list(marker_levels = c("m1"), marker_weights = 1),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaCorr = NULL,
    config = list(
      transforms = list(),
      transforms_spec = list(
        cv_total = list(type = "functional", expr = ~ expit(x))
      ),
      association_plot_payload = list(
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
    call = quote(joinme::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  p <- plot(fit, type = "association", association_term = "cv_total", association_metric = "hazard")
  median_col <- joinme:::.quantile_name_from_prob(0.5)

  expect_equal(range(p$data$x), c(-3, 3), tolerance = 1e-8)
  expect_equal(p$data[[median_col]], 2 * (1 / (1 + exp(-p$data$x))), tolerance = 1e-8)
})
