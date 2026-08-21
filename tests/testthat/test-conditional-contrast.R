.conditional_contrast_test_fit <- function() {
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
  ) # deterministic posterior sample permitting exact paired-contrast assertions

  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "conditional_contrast_mock_fit"),
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
      time_var = "time",
      include_survival = 1L
    ),
    formulaLong = y ~ time + x + trt + (1 | id) + (1 | marker),
    formulaEvent = survival::Surv(event_time, event) ~ x + trt,
    formulaVCov = NULL,
    config = list(tmax = 2, include_survival = TRUE),
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
  ) # fitted-object shell carrying the same design metadata as a sampled model

  list(fit = fit, draws = posterior_draws)
}

test_that("conditional_contrast forms longitudinal differences within posterior draws", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  conditions <- data.frame(
    time = c(0, 1),
    x = c(-0.5, 0.5),
    trt = factor(c("A", "A"), levels = c("A", "B")),
    cond__ = c("early", "late")
  ) # common covariate profiles; trt is intentionally overridden by both groups
  contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = conditions,
    process = "longitudinal",
    method = "posterior_linpred",
    robust = TRUE,
    plot = FALSE
  )

  expected_draws <- fixture$draws[, "beta[4]"] # paired B-minus-A predictor difference in every MCMC draw
  longitudinal <- contrast$longitudinal
  expect_s3_class(contrast, "JoiNMeConditionalContrasts")
  expect_named(contrast, "longitudinal")
  expect_equal(nrow(longitudinal), 2L)
  expect_equal(longitudinal$estimate__, rep(stats::median(expected_draws), 2L))
  expect_equal(longitudinal$se__, rep(stats::sd(expected_draws), 2L))
  expect_equal(longitudinal$lower__, rep(stats::quantile(expected_draws, 0.025, names = FALSE), 2L))
  expect_equal(longitudinal$upper__, rep(stats::quantile(expected_draws, 0.975, names = FALSE), 2L))
  expect_identical(unique(longitudinal$marker__), "Marginal over markers")
  expect_false("trt" %in% names(longitudinal))
  expect_identical(unique(longitudinal$contrast__), "groupA - groupB")
  expect_identical(unique(longitudinal$method__), "posterior_linpred")
})

test_that("conditional_contrast preserves logical predictors from matrix-like fitted data", {
  fitted_matrix <- matrix(
    c(TRUE, FALSE),
    ncol = 1L,
    dimnames = list(NULL, "arm")
  ) # matrix-like stored analysis data that previously lost its name after `$<-`
  process_spec <- list(
    variables = "arm",
    data = fitted_matrix
  ) # minimal process description used by the counterfactual profile builder
  condition_rows <- data.frame(
    cond__ = factor("1", levels = "1")
  ) # one default common condition

  recovered_spec <- .conditional_effect_process_spec(
    formula = y ~ arm,
    data = cbind(y = c(1, 2), fitted_matrix),
    process = "longitudinal"
  ) # public-method preparation must normalise the stored matrix before subsetting

  profile_a <- .conditional_contrast_evaluation_data(
    process_spec,
    condition_rows,
    list(arm = TRUE)
  ) # logical group-A profile
  profile_b <- .conditional_contrast_evaluation_data(
    process_spec,
    condition_rows,
    list(arm = FALSE)
  ) # logical group-B profile

  expect_s3_class(recovered_spec$data, "data.frame")
  expect_identical(recovered_spec$variables, "arm")
  expect_s3_class(profile_a, "data.frame")
  expect_named(profile_a, c("arm", "cond__", "condition_index__"))
  expect_identical(profile_a$arm, TRUE)
  expect_identical(profile_b$arm, FALSE)
  expect_silent(stats::model.matrix(~ arm, data = profile_a))
  expect_silent(stats::model.matrix(~ arm, data = profile_b))
})

test_that("conditional_contrast returns paired MCMC samples when summary is false", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  conditions <- data.frame(time = c(0, 1), cond__ = c("early", "late"))
  raw <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = conditions,
    process = c("longitudinal", "event"),
    method = "posterior_linpred",
    event_scale = "hazard_ratio",
    summary = FALSE
  )

  expect_s3_class(raw, "JoiNMeConditionalContrasts")
  expect_false(attr(raw, "summary"))
  expect_s3_class(raw$longitudinal, "JoiNMeConditionalMCMCSamples")
  expect_s3_class(raw$event, "JoiNMeConditionalMCMCSamples")
  expect_s3_class(raw$longitudinal, "data.frame")
  expect_true(all(c(".draw", ".value", "estimand__") %in% names(raw$longitudinal)))
  expect_equal(nrow(raw$longitudinal), 8L)
  expect_equal(length(unique(raw$longitudinal$estimand__)), 2L)
  expect_equal(
    raw$longitudinal$.value,
    rep(fixture$draws[, "beta[4]"], 2L)
  )
  expect_equal(
    raw$event$.value,
    rep(exp(fixture$draws[, "gamma_w[1,2]"]), 2L)
  )
  expect_false(any(c("estimate__", "se__", "lower__", "upper__") %in% names(raw$event)))
  expect_message(print(raw), "draw-level MCMC samples")
  expect_error(plot(raw), "cannot be plotted directly")
})

test_that("conditional_contrast samples preserve marker marginalisation and event causes", {
  fixture <- .conditional_contrast_test_fit()
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

  marginal <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(time = c(0, 1)),
    process = "longitudinal",
    method = "posterior_linpred",
    longitudinal_estimand = "marginal_marker",
    summary = FALSE
  )$longitudinal
  expect_equal(nrow(marginal), 8L)
  expect_equal(length(unique(marginal$estimand__)), 2L)
  expect_identical(unique(marginal$marker__), "Marginal over markers")

  event <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(x = c(0, 1)),
    process = "event",
    event_scale = "log_hazard_ratio",
    summary = FALSE
  )$event
  expect_identical(unique(event$event_type__), 1:2)
  expect_equal(event$.value[event$event_type__ == 1L], rep(fixture$draws[, "gamma_w[1,2]"], 2L))
  expect_equal(event$.value[event$event_type__ == 2L], rep(fixture$draws[, "gamma_w[2,2]"], 2L))
})

test_that("deterministic event contrasts retain one value per posterior draw", {
  fixture <- .conditional_contrast_test_fit()
  fixture$fit$formulaEvent <- survival::Surv(event_time, event) ~ 1
  fixture$fit$stan_data$p_w <- 0L
  fixture$fit$stan_data$w_cols <- character(0)
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  raw <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(time = c(0, 1)),
    process = c("longitudinal", "event"),
    event_scale = "hazard_ratio",
    summary = FALSE
  )$event

  expect_equal(nrow(raw), 8L)
  expect_equal(length(unique(raw$estimand__)), 2L)
  expect_equal(raw$.value, rep(1, 8L))
})

test_that("conditional_contrast accepts a common fitted subject", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(id = 1, time = c(0, 1), cond__ = c("early", "late")),
    process = "longitudinal",
    method = "posterior_linpred",
    reuse_fitted_re = TRUE,
    plot = FALSE
  )
  expect_s3_class(contrast, "JoiNMeConditionalContrasts")
  expect_true(isTRUE(attr(contrast, "reuse_fitted_re")))
  expect_true(all(contrast$longitudinal$id == 1))
  expect_error(
    conditional_contrast(
      fixture$fit,
      groupA = c(id = 1, trt = "B"),
      groupB = c(id = 2, trt = "A"),
      conditions = data.frame(time = 0),
      process = "longitudinal",
      reuse_fitted_re = TRUE,
      plot = FALSE
    ),
    "identifier cannot differ"
  )
})

test_that("conditional_contrast marginalises markers before longitudinal summary", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  contrast <- conditional_contrast(
    fixture$fit,
    groupA = list(trt = "B", x = 1),
    groupB = list(trt = "A", x = 1),
    conditions = data.frame(time = c(0, 1), cond__ = c("zero", "one")),
    process = "longitudinal",
    method = "posterior_linpred",
    longitudinal_estimand = "marginal_marker",
    plot = FALSE
  )$longitudinal

  expect_equal(nrow(contrast), 2L)
  expect_identical(unique(contrast$marker__), "Marginal over markers")
  expect_equal(contrast$estimate__, rep(stats::median(fixture$draws[, "beta[4]"]), 2L))
})

test_that("conditional_contrast reports paired event hazard ratios and log hazard ratios", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  conditions <- data.frame(x = c(0, 1), cond__ = c("x=0", "x=1")) # common event-regression profiles
  log_contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = conditions,
    process = "survival",
    event_scale = "log_hazard_ratio",
    plot = FALSE
  )$event
  hazard_ratio <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = conditions,
    process = "event",
    event_scale = "hazard_ratio",
    plot = FALSE
  )$event

  expected_log_draws <- fixture$draws[, "gamma_w[1,2]"] # direct event eta_B minus eta_A in each draw
  expect_equal(log_contrast$estimate__, rep(stats::median(expected_log_draws), 2L))
  expect_equal(log_contrast$se__, rep(stats::sd(expected_log_draws), 2L))
  expect_equal(hazard_ratio$estimate__, rep(stats::median(exp(expected_log_draws)), 2L))
  expect_equal(hazard_ratio$se__, rep(stats::sd(exp(expected_log_draws)), 2L))
  expect_identical(unique(hazard_ratio$contrast__), "groupA / groupB")
  expect_identical(unique(hazard_ratio$event_scale__), "hazard_ratio")
})

test_that("conditional_contrast prints and plots process-specific results", {
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(time = c(0, 1), x = 0, cond__ = c("start", "end")),
    process = c("longitudinal", "event"),
    method = "posterior_linpred",
    plot = FALSE
  )
  plots <- plot(contrast, plot = FALSE, condition_variable = "time")

  expect_named(plots, c("longitudinal", "event"))
  expect_s3_class(plots$longitudinal, "ggplot")
  expect_s3_class(plots$event, "ggplot")
  expect_message(returned <- print(contrast), "Comparison:.*trt=B.*trt=A")
  expect_message(print(contrast), "Longitudinal process: 2 conditions")
  expect_identical(returned, contrast)
})

test_that("conditional_contrast plots share the conditional arrangement interface", {
  skip_if_not_installed("patchwork")
  fixture <- .conditional_contrast_test_fit()
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    conditions = data.frame(time = c(0, 1), x = 0, cond__ = c("start", "end")),
    process = c("longitudinal", "event"),
    method = "posterior_linpred",
    plot = FALSE
  )

  arranged <- plot(
    contrast,
    plot = FALSE,
    condition_variable = "time",
    arrange = "grid",
    nrow = 1,
    widths = c(2, 1),
    guides = "collect"
  )

  expect_s3_class(arranged, "patchwork")
  expect_error(
    plot(contrast, plot = FALSE, arrange = "row", ncol = 2),
    "determines the dimensions"
  )
  expect_error(
    plot(contrast, plot = FALSE, arrange = "grid", design = "AB"),
    "arrange = 'design'"
  )
})

test_that("conditional_contrast validates profile names, values, and processes", {
  fixture <- .conditional_contrast_test_fit()

  expect_error(
    conditional_contrast(fixture$fit, groupA = "B", groupB = c(trt = "A"), plot = FALSE),
    "must have a predictor name"
  )
  expect_error(
    conditional_contrast(fixture$fit, groupA = c(unknown = 1), groupB = c(unknown = 0), plot = FALSE),
    "Unknown or unavailable"
  )
  expect_error(
    conditional_contrast(fixture$fit, groupA = c(trt = "A"), groupB = c(trt = "A"), plot = FALSE),
    "must describe different"
  )
  expect_error(
    conditional_contrast(fixture$fit, groupA = c(trt = "B"), groupB = c(trt = "A"), process = "not-a-process", plot = FALSE),
    "Unknown.*process"
  )
})

test_that("conditional_contrast adapts its default to a longitudinal-only fit", {
  fixture <- .conditional_contrast_test_fit()
  fixture$fit$config$include_survival <- FALSE
  fixture$fit$stan_data$include_survival <- 0L
  fixture$fit$formulaEvent <- NULL
  fixture$fit$dataEvent <- NULL
  testthat::local_mocked_bindings(
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      keep <- intersect(variables, colnames(fixture$draws))
      fixture$draws[, keep, drop = FALSE]
    },
    .package = "joinme"
  )

  contrast <- conditional_contrast(
    fixture$fit,
    groupA = c(trt = "B"),
    groupB = c(trt = "A"),
    method = "posterior_linpred",
    plot = FALSE
  )

  expect_named(contrast, "longitudinal")
  expect_error(
    conditional_contrast(
      fixture$fit,
      groupA = c(trt = "B"),
      groupB = c(trt = "A"),
      process = "event",
      plot = FALSE
    ),
    "unavailable for a longitudinal-only"
  )
})
