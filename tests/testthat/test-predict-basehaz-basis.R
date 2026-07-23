test_that("dynamic prediction reuses a fitted natural-spline baseline basis", {
  event_time <- c(0.2, 0.4, 0.6, 0.8, 1)
  fitted_basis <- splines::ns(
    event_time,
    knots = c(0.35, 0.7),
    Boundary.knots = c(0, 1),
    intercept = TRUE
  )
  centre <- seq(0.05, 0.05 * ncol(fitted_basis), by = 0.05)
  object <- list(
    stan_data = list(
      S_event = event_time,
      Kbs = ncol(fitted_basis),
      Bs_event_c = sweep(as.matrix(fitted_basis), 2, centre, "-"),
      Bs_obj = fitted_basis,
      basehaz = "ns",
      basehaz_degree = 3L,
      tmax = 3
    ),
    config = list(),
    call = quote(joinme(id_var = "id", time_var = "time")),
    dataEvent = data.frame(id = seq_along(event_time), time = event_time * 3)
  )

  metadata <- joinme:::.recover_metadata(object, tmax_arg = NULL)
  prediction_time <- c(0.1, 0.3, 0.9)
  prediction_basis <- joinme:::.evaluate_fitted_basehaz_basis(
    object,
    scaled_time = prediction_time,
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
    col_means = numeric(2)
  )

  expect_error(
    joinme:::.evaluate_fitted_basehaz_basis(
      object,
      scaled_time = c(0.2, 0.8),
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

  metadata <- joinme:::.recover_metadata(object, tmax_arg = NULL)
  basis <- joinme:::.evaluate_fitted_basehaz_basis(
    object,
    scaled_time = c(0.2, 0.8),
    basehaz = metadata$basehaz,
    data_event = object$dataEvent[1, , drop = FALSE]
  )

  expect_equal(metadata$basehaz$col_means, c(0, 0.5))
  expect_equal(dim(basis), c(2L, 2L))
  expect_equal(as.numeric(basis), as.numeric(cbind(1, c(0.2, 0.8))))
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
  expect_equal(formula$basehaz_formula, ~ 1 + time, ignore_environment = TRUE)
  expect_length(formula$basehaz_col_means, formula$Kbs)
  expect_equal(formula$basehaz_cols, c("(Intercept)", "time"))

  fit_like <- list(
    stan_data = natural,
    config = list(),
    call = quote(joinme(id_var = "id", time_var = "time")),
    dataEvent = sim$dataEvent
  )
  metadata <- joinme:::.recover_metadata(fit_like, tmax_arg = NULL)
  evaluated <- joinme:::.evaluate_fitted_basehaz_basis(
    fit_like,
    scaled_time = c(0.1, 0.5, 0.9),
    basehaz = metadata$basehaz,
    data_event = sim$dataEvent[1, , drop = FALSE]
  )
  expect_equal(dim(evaluated), c(3L, natural$Kbs))
  expect_equal(metadata$basehaz$col_means, natural$basehaz_col_means)
})
