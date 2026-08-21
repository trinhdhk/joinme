test_that("dynamic prediction reuses a fitted natural-spline baseline basis", {
  scaled_event_time <- c(0.2, 0.4, 0.6, 0.8, 1)
  event_time <- scaled_event_time * 3
  fitted_basis <- splines::ns(
    event_time,
    knots = c(1.05, 2.1),
    Boundary.knots = c(0, 3),
    intercept = TRUE
  )
  centre <- seq(0.05, 0.05 * ncol(fitted_basis), by = 0.05)
  object <- list(
    stan_data = list(
      S_event = scaled_event_time,
      Kbs = ncol(fitted_basis),
      Bs_event_c = sweep(as.matrix(fitted_basis), 2, centre, "-"),
      Bs_obj = fitted_basis,
      basehaz = "ns",
      basehaz_degree = 3L,
      tmax = 3
    ),
    config = list(),
    call = quote(joinme(id_var = "id", time_var = "time")),
    dataEvent = data.frame(id = seq_along(event_time), time = event_time)
  )

  metadata <- .recover_metadata(object, tmax_arg = NULL)
  prediction_time <- c(0.3, 0.9, 2.7)
  prediction_basis <- .evaluate_fitted_basehaz_basis(
    object,
    original_time = prediction_time,
    basehaz = metadata$basehaz,
    data_event = object$dataEvent[1, , drop = FALSE]
  )

  expect_identical(metadata$basehaz$type, "ns")
  expect_identical(metadata$basehaz$basis_object, fitted_basis)
  expect_equal(metadata$basehaz$col_means, centre)
  expect_equal(dim(prediction_basis), c(3L, ncol(fitted_basis)))
  expect_equal(
    prediction_basis,
    as.matrix(predict(fitted_basis, newx = prediction_time))
  )
})

test_that("baseline-basis validation reports a fitted-column mismatch early", {
  object <- list(
    stan_data = list(Kbs = 2L),
    call = quote(joinme(id_var = "id", time_var = "time"))
  )
  metadata <- list(
    type = "bs",
    knots = c(0.25, 0.75),
    degree = 3L,
    formula = NULL,
    basis_object = NULL,
    time_scale = 1,
    col_means = numeric(2)
  )

  expect_error(
    .evaluate_fitted_basehaz_basis(
      object,
      original_time = c(0.2, 0.8),
      basehaz = metadata
    ),
    "incompatible baseline-hazard basis"
  )
})

test_that("formula baseline metadata retains its fitted centring constants", {
  object <- list(
    stan_data = list(
      S_event = c(0.25, 0.5, 1),
      Kbs = 2L,
      Bs_event_c = cbind(1, c(-0.25, 0, 0.5)),
      basehaz = "formula",
      basehaz_formula = ~ 1 + time,
      basehaz_col_means = c(0, 0.5),
      tmax = 4
    ),
    config = list(),
    call = quote(joinme(id_var = "id", time_var = "time")),
    dataEvent = data.frame(id = 1:3, time = c(1, 2, 4))
  )

  metadata <- .recover_metadata(object, tmax_arg = NULL)
  basis <- .evaluate_fitted_basehaz_basis(
    object,
    original_time = c(0.8, 3.2),
    basehaz = metadata$basehaz,
    data_event = object$dataEvent[1, , drop = FALSE]
  )

  expect_equal(metadata$basehaz$col_means, c(0, 0.5))
  expect_equal(dim(basis), c(2L, 2L))
  expect_equal(as.numeric(basis), as.numeric(cbind(1, c(0.8, 3.2))))
})

test_that("formula baselines recognise a Surv clock distinct from longitudinal time", {
  training_event <- data.frame(
    id = 1:4,
    visit_time = c(1, 2, 3, 4),
    time = 0
  )
  baseline_template <- .make_model_matrix_template(
    ~ 1 + splines::ns(visit_time, knots = 2, Boundary.knots = c(0, 4)),
    training_event
  )
  object <- list(
    stan_data = list(
      Kbs = length(baseline_template$columns),
      event_time_vars = c("time", "visit_time")
    ),
    call = quote(joinme(id_var = "id", time_var = "time"))
  )
  metadata <- list(
    type = "formula",
    formula = baseline_template,
    basis_object = NULL,
    time_scale = 4
  )
  prediction_time <- c(0.5, 2.5, 3.5)

  evaluated <- .evaluate_fitted_basehaz_basis(
    object,
    original_time = prediction_time,
    basehaz = metadata,
    data_event = training_event[1L, , drop = FALSE]
  )
  expected_data <- training_event[rep(1L, length(prediction_time)), , drop = FALSE]
  expected_data$time <- prediction_time
  expected_data$visit_time <- prediction_time

  expect_equal(evaluated, .mm(baseline_template, expected_data))
})

test_that("standata retains baseline metadata required by dynamic prediction", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 3,
    D = 2,
    n_t = 2,
    seed = 2207,
    include_marker_only = TRUE
  )
  common <- list(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent
  )

  natural <- do.call(joinme_standata, c(common, list(
    basehaz = "ns",
    basehaz_n_knots = 2L
  )))
  formula <- do.call(joinme_standata, c(common, list(
    basehaz = "formula",
    basehaz_formula = ~ 1 + time
  )))

  expect_s3_class(natural$Bs_obj, "ns")
  expect_length(natural$basehaz_col_means, natural$Kbs)
  expect_true(all(is.finite(natural$basehaz_col_means)))
  expect_equal(natural$basehaz_cols, paste0("basis_", seq_len(natural$Kbs)))
  expect_s3_class(formula$basehaz_formula, "JoiNMe_mm")
  expect_equal(formula$basehaz_formula$formula, ~ 1 + time, ignore_environment = TRUE)
  expect_length(formula$basehaz_col_means, formula$Kbs)
  expect_equal(formula$basehaz_cols, c("(Intercept)", "time"))

  fit_like <- list(
    stan_data = natural,
    config = list(),
    call = quote(joinme(id_var = "id", time_var = "time")),
    dataEvent = sim$dataEvent
  )
  metadata <- .recover_metadata(fit_like, tmax_arg = NULL)
  evaluated <- .evaluate_fitted_basehaz_basis(
    fit_like,
    original_time = c(0.1, 0.5, 0.9) * natural$tmax,
    basehaz = metadata$basehaz,
    data_event = sim$dataEvent[1, , drop = FALSE]
  )
  expect_equal(dim(evaluated), c(3L, natural$Kbs))
  expect_equal(metadata$basehaz$col_means, natural$basehaz_col_means)
})

test_that("baseline-hazard coefficient labels follow their statistical representation", {
  expect_equal(
    .basehaz_term_labels("bs", 3L, c("1", "2", "3")),
    c("basis_1", "basis_2", "basis_3")
  )
  expect_equal(
    .basehaz_term_labels("ns", 2L, c("1", "2")),
    c("basis_1", "basis_2")
  )
  expect_equal(
    .basehaz_term_labels("formula", 3L, c("(Intercept)", "time", "I(time^2)")),
    c("(Intercept)", "time", "I(time^2)")
  )
  expect_equal(
    .basehaz_term_labels("pwlin", 3L, c("1", "2", "3")),
    c("segment_1", "segment_2", "segment_3")
  )
  expect_equal(
    .basehaz_term_labels("pwlin", 2L, c("before_2", "after_2")),
    c("before_2", "after_2")
  )
})
