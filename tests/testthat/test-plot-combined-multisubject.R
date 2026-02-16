testthat::test_that("combined multi-subject returns per-subject combined plots", {
  testthat::local_mocked_bindings(
    .combine_plot_grid = function(plots, ncol = NULL, fallback = c("flat", "input")) {
      ggplot2::ggplot()
    },
    .package = "joinme"
  )

  ids <- c(1, 2)
  time_grid <- seq(0, 1, length.out = 5)

  quant_long <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = time_grid,
      marker = "m1",
      q2.5 = seq(0.1, 0.5, length.out = length(time_grid)) + 0.05 * id,
      q50 = seq(0.2, 0.6, length.out = length(time_grid)) + 0.05 * id,
      q97.5 = seq(0.3, 0.7, length.out = length(time_grid)) + 0.05 * id,
      mean = seq(0.2, 0.6, length.out = length(time_grid)) + 0.05 * id,
      sd = rep(0.05, length(time_grid))
    )
  }))

  quant_surv <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = time_grid,
      q2.5 = pmax(0, 0.8 - time_grid - 0.03 * id),
      q50 = pmax(0, 0.9 - time_grid - 0.02 * id),
      q97.5 = pmax(0, 1.0 - time_grid - 0.01 * id),
      mean = pmax(0, 0.9 - time_grid - 0.02 * id),
      sd = rep(0.03, length(time_grid))
    )
  }))

  data_long <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = seq(0, 1, length.out = 3),
      marker = "m1",
      y = seq(0.2, 0.4, length.out = 3)
    )
  }))

  pred <- joinme::JoinMeDynPred$new(
    predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
    quantiles = list(
      longitudinal = quant_long,
      longitudinal_fitted = NULL,
      longitudinal_marker_pop = NULL,
      longitudinal_overall_pop = NULL,
      survival = quant_surv,
      cumhaz = NULL
    ),
    draws = list(longitudinal = NULL, longitudinal_fitted = NULL, survival = list(), cumhaz = list()),
    data = list(longitudinal = data_long, event = data.frame(id = ids, time = 1, event = 0L)),
    metadata = list(
      scale = "epred",
      ci_levels = 0.95,
      conditioning_time = 0.5,
      conditioning_time_by_id = c(`1` = 0.5, `2` = 0.5),
      id_var = "id",
      time_var = "time",
      marker_var = "marker",
      response_var = "y"
    ),
    call = quote(predict(fit_obj)),
    tmax = 1,
    n_samples = 20
  )

  out <- plot(
    pred,
    which = c("longitudinal", "survival"),
    subject = ids,
    combined = TRUE,
    ci_levels = 0.95,
    smooth_trajectory = FALSE,
    show_data = FALSE,
    facet_by = "none"
  )

  testthat::expect_type(out, "list")
  testthat::expect_setequal(names(out), as.character(ids))
  testthat::expect_true(all(vapply(out, function(x) inherits(x, "gg"), logical(1))))
})

testthat::test_that("combined multi-subject fallback preserves per-subject structure", {
  testthat::local_mocked_bindings(
    .combine_plot_grid = function(plots, ncol = NULL, fallback = c("flat", "input")) {
      fallback <- match.arg(fallback)
      if (identical(fallback, "input")) return(plots)
      list(flat = ggplot2::ggplot())
    },
    .package = "joinme"
  )

  ids <- c(1, 2)
  time_grid <- seq(0, 1, length.out = 5)

  quant_long <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = time_grid,
      marker = "m1",
      q2.5 = seq(0.1, 0.5, length.out = length(time_grid)) + 0.05 * id,
      q50 = seq(0.2, 0.6, length.out = length(time_grid)) + 0.05 * id,
      q97.5 = seq(0.3, 0.7, length.out = length(time_grid)) + 0.05 * id,
      mean = seq(0.2, 0.6, length.out = length(time_grid)) + 0.05 * id,
      sd = rep(0.05, length(time_grid))
    )
  }))

  quant_surv <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = time_grid,
      q2.5 = pmax(0, 0.8 - time_grid - 0.03 * id),
      q50 = pmax(0, 0.9 - time_grid - 0.02 * id),
      q97.5 = pmax(0, 1.0 - time_grid - 0.01 * id),
      mean = pmax(0, 0.9 - time_grid - 0.02 * id),
      sd = rep(0.03, length(time_grid))
    )
  }))

  data_long <- do.call(rbind, lapply(ids, function(id) {
    data.frame(
      id = id,
      time = seq(0, 1, length.out = 3),
      marker = "m1",
      y = seq(0.2, 0.4, length.out = 3)
    )
  }))

  pred <- joinme::JoinMeDynPred$new(
    predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
    quantiles = list(
      longitudinal = quant_long,
      longitudinal_fitted = NULL,
      longitudinal_marker_pop = NULL,
      longitudinal_overall_pop = NULL,
      survival = quant_surv,
      cumhaz = NULL
    ),
    draws = list(longitudinal = NULL, longitudinal_fitted = NULL, survival = list(), cumhaz = list()),
    data = list(longitudinal = data_long, event = data.frame(id = ids, time = 1, event = 0L)),
    metadata = list(
      scale = "epred",
      ci_levels = 0.95,
      conditioning_time = 0.5,
      conditioning_time_by_id = c(`1` = 0.5, `2` = 0.5),
      id_var = "id",
      time_var = "time",
      marker_var = "marker",
      response_var = "y"
    ),
    call = quote(predict(fit_obj)),
    tmax = 1,
    n_samples = 20
  )

  out <- plot(
    pred,
    which = c("longitudinal", "survival"),
    subject = ids,
    combined = TRUE,
    ci_levels = 0.95,
    smooth_trajectory = FALSE,
    show_data = FALSE,
    facet_by = "none"
  )

  testthat::expect_type(out, "list")
  testthat::expect_setequal(names(out), as.character(ids))
  testthat::expect_true(all(vapply(out, is.list, logical(1))))
  testthat::expect_true(all(vapply(out, function(x) all(c("longitudinal", "survival") %in% names(x)), logical(1))))
})
