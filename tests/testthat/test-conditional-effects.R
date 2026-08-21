.conditional_effect_test_fit <- function() {
  posterior_draws <- cbind(
    `beta[1]` = c(0.1, 0.2, 0.3, 0.4),
    `beta[2]` = c(0.8, 0.9, 1.0, 1.1),
    `beta[3]` = c(0.2, 0.3, 0.4, 0.5),
    `beta[4]` = c(-0.1, -0.2, -0.3, -0.4),
    `v_marker[1,1]` = rep(-1, 4),
    `v_marker[2,1]` = rep(1, 4),
    `u_id[1,1]` = rep(0.75, 4),
    `u_id[2,1]` = rep(-0.25, 4),
    `gamma_w[1,1]` = c(0.2, 0.3, 0.4, 0.5),
    `gamma_w[1,2]` = c(-0.4, -0.3, -0.2, -0.1)
  )

  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "conditional_effect_mock_fit"),
    stan_data = list(
      P = 4L,
      D = 2L,
      R_id = 1L,
      R_mk = 1L,
      Q_idm = 0L,
      n_id = 2L,
      p_w = 2L,
      K_event = 1L,
      w_cols = c("x", "trtB"),
      marker_levels = c("m1", "m2"),
      link_long = c(3L, 3L),
      inv_link_n_ops = c(0L, 0L),
      inv_link_n_const = c(0L, 0L),
      inv_link_ops = matrix(0L, nrow = 2L, ncol = 1L),
      inv_link_const = matrix(0, nrow = 2L, ncol = 1L),
      time_var = "time"
    ),
    formulaLong = y ~ time + x + trt + (1 | id) + (1 | marker),
    formulaEvent = survival::Surv(event_time, event) ~ x + trt,
    formulaVCov = NULL,
    config = list(tmax = 2),
    call = quote(joinme(
      formulaLong = y ~ time + x + trt + (1 | id) + (1 | marker),
      formulaEvent = survival::Surv(event_time, event) ~ x + trt,
      id_var = "id",
      time_var = "time",
      marker_var = "marker"
    )),
    tmax = 2,
    dataLong = data.frame(
      id = rep(1:2, each = 4),
      time = rep(c(0, 1), 4),
      marker = factor(rep(c("m1", "m2"), each = 2, times = 2), levels = c("m1", "m2")),
      x = rep(c(-1, 1), each = 4),
      trt = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
      y = seq_len(8) / 10
    ),
    dataEvent = data.frame(
      id = 1:4,
      event_time = c(1, 1.5, 2, 2),
      event = c(0L, 1L, 0L, 1L),
      x = c(-1, 0, 1, 2),
      trt = factor(c("A", "A", "B", "B"), levels = c("A", "B"))
    )
  )

  list(fit = fit, draws = posterior_draws)
}

test_that("conditional_effects can reuse paired fitted subject effects", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  fitted_subject_effects <- conditional_effects(
    fixture$fit,
    effects = "time",
    conditions = data.frame(id = c(1, 2), cond__ = c("subject 1", "subject 2")),
    process = "longitudinal",
    method = "posterior_linpred",
    longitudinal_estimand = "population",
    markers = "m1",
    resolution = 2,
    reuse_fitted_re = TRUE,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time

  estimates_by_subject <- split(fitted_subject_effects$estimate__, fitted_subject_effects$id)
  expect_equal(estimates_by_subject[["1"]] - estimates_by_subject[["2"]], rep(1, 2))
  expect_true(isTRUE(attr(
    conditional_effects(
      fixture$fit,
      effects = "time",
      conditions = data.frame(id = 1, cond__ = "subject 1"),
      process = "longitudinal",
      method = "posterior_linpred",
      markers = "m1",
      resolution = 2,
      reuse_fitted_re = TRUE,
      plot = FALSE
    ),
    "reuse_fitted_re"
  )))
  expect_error(
    conditional_effects(
      fixture$fit,
      effects = "time",
      conditions = data.frame(id = 99),
      process = "longitudinal",
      reuse_fitted_re = TRUE,
      plot = FALSE
    ),
    "Unknown id value"
  )
})

test_that("conditional_effects returns both joint-model processes as brms-compatible data", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  conditions <- make_conditions(fixture$fit$dataEvent, vars = "trt")
  effects <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    conditions = conditions,
    process = c("longitudinal", "event"),
    resolution = 5,
    plot = FALSE
  )

  expect_s3_class(effects, "JoiNMeConditionalEffects")
  expect_named(effects, c("longitudinal", "event"))
  expect_s3_class(effects$longitudinal, "brms_conditional_effects")
  expect_s3_class(effects$event, "brms_conditional_effects")

  longitudinal_data <- effects$longitudinal$time
  event_data <- effects$event$x
  required_columns <- c("effect1__", "estimate__", "se__", "lower__", "upper__", "cond__", "process__")
  expect_true(all(required_columns %in% names(longitudinal_data)))
  expect_true(all(required_columns %in% names(event_data)))
  expect_equal(nrow(longitudinal_data), 5L * 2L)
  expect_equal(nrow(event_data), 5L * 2L)
  expect_identical(unique(longitudinal_data$marker__), "Marginal over markers")
  expect_identical(unique(longitudinal_data$longitudinal_estimand__), "population")
  expect_true(all(event_data$estimate__ > 0))
  expect_identical(attr(longitudinal_data, "response"), "Population expected longitudinal response")
  expect_identical(attr(event_data, "response"), "Relative event hazard")

  both_alias <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    process = "both",
    resolution = 2,
    plot = FALSE
  )
  expect_named(both_alias, c("longitudinal", "event"))
})

test_that("conditional_effects returns draw-level MCMC samples when summary is false", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  raw <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    process = c("longitudinal", "event"),
    markers = "m1",
    resolution = 2,
    summary = FALSE
  )
  summarised <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    process = c("longitudinal", "event"),
    markers = "m1",
    resolution = 2,
    plot = FALSE
  )

  expect_s3_class(raw, "JoiNMeConditionalEffects")
  expect_false(attr(raw, "summary"))
  expect_s3_class(raw$longitudinal$time, "JoiNMeConditionalMCMCSamples")
  expect_s3_class(raw$event$x, "JoiNMeConditionalMCMCSamples")
  expect_s3_class(raw$longitudinal$time, "data.frame")
  expect_true(all(c(".draw", ".value", "estimand__") %in% names(raw$longitudinal$time)))
  expect_equal(sort(unique(raw$longitudinal$time$.draw)), 1:4)
  expect_equal(length(unique(raw$longitudinal$time$estimand__)), nrow(summarised$longitudinal$time))
  expect_false(any(c("estimate__", "se__", "lower__", "upper__") %in% names(raw$longitudinal$time)))
  longitudinal_medians <- vapply(unique(raw$longitudinal$time$estimand__), function(estimand) {
    stats::median(raw$longitudinal$time$.value[raw$longitudinal$time$estimand__ == estimand])
  }, numeric(1)) # posterior centre for every tidy longitudinal estimand
  event_medians <- vapply(unique(raw$event$x$estimand__), function(estimand) {
    stats::median(raw$event$x$.value[raw$event$x$estimand__ == estimand])
  }, numeric(1)) # posterior centre for every tidy event estimand
  expect_equal(
    unname(longitudinal_medians),
    summarised$longitudinal$time$estimate__
  )
  expect_equal(
    unname(event_medians),
    summarised$event$x$estimate__
  )
  expect_message(print(raw), "draw-level MCMC samples")
  expect_error(plot(raw), "cannot be plotted directly")
  expect_error(
    conditional_effects(fixture$fit, effects = "time", process = "longitudinal", summary = NA),
    "summary.*TRUE or FALSE"
  )
})

test_that("conditional_effect draw samples retain fitted-subject and competing-event labels", {
  fixture <- .conditional_effect_test_fit()
  fixture$fit$stan_data$K_event <- 2L
  fixture$draws <- cbind(
    fixture$draws,
    `gamma_w[2,1]` = c(-0.2, -0.3, -0.4, -0.5),
    `gamma_w[2,2]` = c(0.4, 0.3, 0.2, 0.1)
  )
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  subject_draws <- conditional_effects(
    fixture$fit,
    effects = "time",
    conditions = data.frame(id = c(1, 2), cond__ = c("subject 1", "subject 2")),
    process = "longitudinal",
    method = "posterior_linpred",
    markers = "m1",
    resolution = 2,
    reuse_fitted_re = TRUE,
    summary = FALSE
  )$longitudinal$time
  expect_true(all(c(1, 2) %in% subject_draws$id))
  subject_one <- subject_draws[subject_draws$id == 1, , drop = FALSE] # tidy draws for the first fitted subject
  subject_two <- subject_draws[subject_draws$id == 2, , drop = FALSE] # corresponding draws for the second fitted subject
  expect_equal(
    subject_one$.value - subject_two$.value,
    rep(1, nrow(subject_one))
  )

  event_draws <- conditional_effects(
    fixture$fit,
    effects = "x",
    process = "event",
    resolution = 2,
    summary = FALSE
  )$event$x
  expect_identical(unique(event_draws$event_type__), 1:2)
  expect_equal(length(unique(event_draws$estimand__)), 4L)
  expect_equal(nrow(event_draws), 16L)
})

test_that("conditional-effects data print with cli bullets", {
  effects <- structure(
    list(
      longitudinal = structure(
        list(arm = data.frame(effect1__ = c(FALSE, TRUE))),
        class = c("brms_conditional_effects", "list")
      )
    ),
    class = c("JoiNMeConditionalEffects", "list"),
    longitudinal_estimand = "population"
  )

  expect_message(
    returned <- print(effects),
    "Longitudinal estimand: population"
  )
  expect_message(print(effects), "Longitudinal process: 1 effect")
  expect_identical(returned, effects)
})

test_that("conditional_effects can return either process and custom focal values", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  longitudinal <- conditional_effects(
    fixture$fit,
    effects = "time:trt",
    process = "longitudinal",
    int_conditions = list(time = c(0, 1, 2), trt = c("A", "B")),
    markers = "m2",
    method = "posterior_linpred",
    plot = FALSE
  )
  expect_named(longitudinal, "longitudinal")
  expect_identical(unique(longitudinal$longitudinal[[1]]$marker__), "Marginal over markers")
  expect_equal(nrow(longitudinal$longitudinal[[1]]), 6L)
  expect_identical(attr(longitudinal$longitudinal[[1]], "response"), "Population longitudinal linear predictor")

  event <- conditional_effects(
    fixture$fit,
    effects = "trt",
    process = "event",
    event_scale = "log_hazard_ratio",
    plot = FALSE
  )
  expect_named(event, "event")
  expect_identical(attr(event$event$trt, "response"), "Log relative event hazard")
  expect_equal(sort(unique(as.character(event$event$trt$effect1__))), c("A", "B"))
})

test_that("longitudinal conditional effects distinguish population, marker, and marker-marginal estimands", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  # Population predictions exclude v_marker and average the marker-specific
  # inverse-link trajectories within every draw, leaving one curve.
  population <- conditional_effects(
    fixture$fit,
    effects = "time",
    process = "longitudinal",
    longitudinal_estimand = "population",
    resolution = 4,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time
  expect_equal(nrow(population), 4L)
  expect_identical(unique(population$marker__), "Marginal over markers")
  expect_identical(attr(population, "longitudinal_estimand"), "population")

  # Marker predictions add the fitted marker-level random effect while still
  # excluding all subject and marker-by-subject random effects.
  marker_linear <- conditional_effects(
    fixture$fit,
    effects = "time",
    process = "longitudinal",
    longitudinal_estimand = "marker",
    method = "posterior_linpred",
    resolution = 4,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time
  marker_linear_by_level <- split(marker_linear$estimate__, marker_linear$marker__)
  expect_equal(marker_linear_by_level$m2 - marker_linear_by_level$m1, rep(2, 4))

  # Marginal expected responses are averaged after marker-specific inverse-link
  # transformation inside each draw. Posterior means consequently equal the
  # arithmetic mean of marker-specific posterior means, not the population
  # expected response obtained by setting the marker deviation to zero.
  marker_response <- conditional_effects(
    fixture$fit,
    effects = "time",
    process = "longitudinal",
    longitudinal_estimand = "marker",
    resolution = 4,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time
  marginal_response <- conditional_effects(
    fixture$fit,
    effects = "time",
    process = "longitudinal",
    longitudinal_estimand = "marginal_marker",
    resolution = 4,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time
  marker_response_by_level <- split(marker_response$estimate__, marker_response$marker__)
  expected_marginal_mean <- rowMeans(do.call(cbind, marker_response_by_level))

  expect_equal(marginal_response$estimate__, expected_marginal_mean)
  expect_false(isTRUE(all.equal(marginal_response$estimate__, population$estimate__)))
  expect_equal(nrow(marginal_response), 4L)
  expect_identical(unique(marginal_response$marker__), "Marginal over markers")
  expect_identical(unique(marginal_response$longitudinal_estimand__), "marginal_marker")
  expect_equal(nrow(attr(marginal_response, "points")), 0L)
  expect_identical(attr(marginal_response, "response"), "Marker-marginal expected longitudinal response")

  marginal_m2 <- conditional_effects(
    fixture$fit,
    effects = "time",
    process = "longitudinal",
    longitudinal_estimand = "marginal_marker",
    markers = "m2",
    resolution = 4,
    robust = FALSE,
    plot = FALSE
  )$longitudinal$time
  expect_equal(marginal_m2$estimate__, marker_response_by_level$m2)
  expect_s3_class(
    plot(
      structure(
        list(longitudinal = structure(list(time = marginal_response), class = c("brms_conditional_effects", "list"))),
        class = c("JoiNMeConditionalEffects", "list")
      ),
      plot = FALSE
    )$longitudinal$time,
    "ggplot"
  )
})

test_that("conditional_effects delegates plotting to the brms data contract", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  effects <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    longitudinal_estimand = "marker",
    resolution = 4,
    plot = FALSE
  )
  plots <- plot(effects, plot = FALSE, points = TRUE)

  expect_named(plots, c("longitudinal", "event"))
  expect_s3_class(plots$longitudinal$time, "ggplot")
  expect_s3_class(plots$event$x, "ggplot")
  expect_true(nrow(attr(effects$longitudinal$time, "points")) > 0L)
  expect_equal(nrow(attr(effects$event$x, "points")), 0L)
})

test_that("conditional_effect plots support grid, row, column, and design arrangements", {
  skip_if_not_installed("patchwork")
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  effects <- conditional_effects(
    fixture$fit,
    effects = list(longitudinal = "time", event = "x"),
    process = c("longitudinal", "event"),
    resolution = 4,
    plot = FALSE
  )

  grid_plot <- plot(effects, plot = FALSE, arrange = "grid", ncol = 2, guides = "collect")
  row_plot <- plot(effects, plot = FALSE, arrange = "row")
  column_plot <- plot(effects, plot = FALSE, arrange = "column")
  design_plot <- plot(effects, plot = FALSE, arrange = "design", design = "AB")

  expect_s3_class(grid_plot, "patchwork")
  expect_s3_class(row_plot, "patchwork")
  expect_s3_class(column_plot, "patchwork")
  expect_s3_class(design_plot, "patchwork")
  expect_error(
    plot(effects, plot = FALSE, ncol = 2),
    "require an arranged conditional display"
  )
  expect_error(
    plot(effects, plot = FALSE, arrange = "design"),
    "requires.*design"
  )
})

test_that("conditional_effects validates process-specific effects and condition labels", {
  fixture <- .conditional_effect_test_fit()

  expect_error(
    conditional_effects(fixture$fit, effects = "not_in_model", process = "event", plot = FALSE),
    "None of the requested"
  )
  expect_error(
    conditional_effects(
      fixture$fit,
      effects = "x",
      conditions = data.frame(x = c(0, 1), cond__ = c("same", "same")),
      process = "event",
      plot = FALSE
    ),
    "unique"
  )
})

test_that("conditional effect grids retain fitted categorical design levels", {
  cast_character <- .conditional_effect_cast("active", c("control", "active", "control"))
  expect_s3_class(cast_character, "factor")
  expect_equal(levels(cast_character), c("active", "control"))

  observed_factor <- factor(c("A", "B", "A"), levels = c("A", "B"))
  contrasts(observed_factor) <- stats::contr.sum(2)
  cast_factor <- .conditional_effect_cast("B", observed_factor)
  expect_equal(levels(cast_factor), levels(observed_factor))
  expect_equal(contrasts(cast_factor), contrasts(observed_factor))
})

test_that("conditional effects support brms-style two-dimensional surfaces", {
  fixture <- .conditional_effect_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  effects <- conditional_effects(
    fixture$fit,
    effects = "time:x",
    process = "longitudinal",
    markers = "m1",
    resolution = 3,
    surface = TRUE,
    plot = FALSE
  )
  surface_data <- effects$longitudinal[[1L]]
  expect_true(isTRUE(attr(surface_data, "surface")))
  expect_equal(nrow(surface_data), 9L)
  expect_s3_class(plot(effects, plot = FALSE, stype = "raster")$longitudinal[[1L]], "ggplot")
})
