test_that("standata evaluates spline time terms on original study time", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 4,
    D = 2,
    n_t = 3,
    seed = 3101,
    include_marker_only = TRUE
  )

  formulaLong <- y ~ 1 + splines::ns(time, df = 3) + x1 +
    (0 + splines::bs(time, df = 3) | id) +
    (0 + x1 + (0 + splines::ns(time, df = 3) | id) | marker)

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total")
  )

  expect_s3_class(sd$design_templates$fixed, "JoiNMe_mm")
  expect_true(sd$P > 0L)
  expect_true(sd$R_id > 0L)
  expect_true(sd$Q_idm > 0L)

  ordered_long <- sim$dataLong[
    order(sim$dataLong$id, sim$dataLong$marker, sim$dataLong$time),
    ,
    drop = FALSE
  ]
  fixed_formula <- reformulas::nobars(reformulas::expandDoubleVerts(formulaLong))
  expected_fixed <- stats::model.matrix(fixed_formula, data = ordered_long)
  expect_equal(unname(sd$X_obs), unname(expected_fixed), tolerance = 1e-12)

  probe <- ordered_long[seq_len(min(6L, nrow(ordered_long))), , drop = FALSE]
  probe[["time"]] <- seq(0, max(sim$dataLong$time), length.out = nrow(probe))
  probe_matrix <- .mm(sd$design_templates$fixed, probe)
  expected_probe <- stats::model.matrix(
    stats::delete.response(stats::terms(stats::model.frame(fixed_formula, data = ordered_long))),
    data = probe
  )
  expect_equal(unname(probe_matrix), unname(expected_probe), tolerance = 1e-12)
})

test_that("explicit spline knots retain their original-time interpretation", {
  sim <- simulate_joinme(
    formulaLong = y ~ splines::ns(time, knots = 1, Boundary.knots = c(0, 5)) +
      (1 + time | id) + (1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    n_id = 6,
    families = rep("gaussian", 2),
    time_cens = 8,
    times_obs = 0:5,
    assoc = "cv_mean",
    truth = jm_truth(assoc_coef = list(slope = c(cv_mean = 0))),
    seed = 31011,
    use_mirai = FALSE
  )
  sim$dataEvent$time <- 8
  sim$dataEvent$event <- 0L

  sd <- joinme_standata(
    formulaLong = sim$truth$formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = sim$truth$formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = "cv_mean",
    families = sim$marker_info$families
  )

  ordered_long <- sim$dataLong[
    order(sim$dataLong$id, sim$dataLong$marker, sim$dataLong$time),
    ,
    drop = FALSE
  ]
  fixed_formula <- reformulas::nobars(
    reformulas::expandDoubleVerts(sim$truth$formulaLong)
  )
  expected_original <- stats::model.matrix(fixed_formula, data = ordered_long)
  scaled_long <- ordered_long
  scaled_long$time <- scaled_long$time / sd$tmax
  expected_scaled <- stats::model.matrix(fixed_formula, data = scaled_long)

  expect_equal(unname(sd$X_obs), unname(expected_original), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(sd$X_obs), unname(expected_scaled))))

  event_design_data <- sim$dataEvent
  event_design_data$time <- sd$S_event * sd$tmax
  expected_event <- .mm(sd$design_templates$fixed, event_design_data)
  expect_equal(unname(sd$X_event_now), unname(expected_event), tolerance = 1e-12)

  event_forward_data <- sim$dataEvent
  event_forward_data$time <- (sd$S_event + sd$eps_fd) * sd$tmax
  expected_event_forward <- .mm(sd$design_templates$fixed, event_forward_data)
  expect_equal(unname(sd$X_event_fwd), unname(expected_event_forward), tolerance = 1e-12)
  expect_true(all(event_forward_data$time > event_design_data$time))

  dynamic_event_ordinates <- c(0.125, 0.5, 0.875)
  dynamic_original_times <- .event_ordinate_to_original_time(
    dynamic_event_ordinates,
    sd$tmax
  )
  dynamic_design_data <- sim$dataEvent[rep(1L, length(dynamic_event_ordinates)), , drop = FALSE]
  dynamic_design_data$time <- dynamic_original_times
  dynamic_design <- .mm(sd$design_templates$fixed, dynamic_design_data)
  wrongly_scaled_data <- dynamic_design_data
  wrongly_scaled_data$time <- dynamic_event_ordinates
  wrongly_scaled_design <- .mm(sd$design_templates$fixed, wrongly_scaled_data)
  expect_false(isTRUE(all.equal(unname(dynamic_design), unname(wrongly_scaled_design))))
  expect_equal(dynamic_original_times, dynamic_event_ordinates * 8)
})

test_that("every formula-derived submodel uses original study time", {
  event_data <- data.frame(
    id = seq_len(4L),
    event_time = c(4, 5, 6, 8),
    event = 0L,
    treatment = c(0, 1, 0, 1)
  )
  longitudinal_data <- do.call(rbind, lapply(seq_len(nrow(event_data)), function(subject) {
    do.call(rbind, lapply(c("m1", "m2"), function(marker_name) {
      data.frame(
        id = subject,
        marker = marker_name,
        time = c(0, 2, 4),
        treatment = event_data$treatment[[subject]],
        y = 0
      )
    }))
  }))

  event_formula <- survival::Surv(event_time, event) ~
    treatment * splines::ns(event_time, knots = 3, Boundary.knots = c(0, 8))
  covariance_formula <- list(
    sd = ~ splines::ns(event_time, knots = 3, Boundary.knots = c(0, 8)),
    corr = ~ treatment:splines::ns(event_time, knots = 3, Boundary.knots = c(0, 8))
  )
  distribution_formula <- list(
    sigma ~ splines::ns(time, knots = 1.5, Boundary.knots = c(0, 8))
  )

  standata <- joinme_standata(
    formulaLong = y ~ treatment + splines::ns(time, knots = 1.5, Boundary.knots = c(0, 8)) +
      (1 + time | id) + (1 + (1 + time | id) | marker),
    dataLong = longitudinal_data,
    formulaEvent = event_formula,
    dataEvent = event_data,
    formulaVCov = covariance_formula,
    formulaDist = distribution_formula,
    families = c("gaussian", "gaussian"),
    assoc = "cv_mean",
    basehaz = "formula",
    basehaz_formula = ~ 1 + splines::ns(event_time, knots = 3, Boundary.knots = c(0, 8)),
    quadrature_nodes = 7L
  )

  expect_identical(standata$event_time_vars, c("time", "event_time"))
  expect_s3_class(standata$design_templates$event, "JoiNMe_mm")
  expect_s3_class(standata$design_templates$vcov$sd, "JoiNMe_mm")
  expect_s3_class(standata$design_templates$distributional$sigma$default, "JoiNMe_mm")

  expected_event <- .mm_event(standata$design_templates$event, event_data)
  expect_equal(unname(standata$W), unname(expected_event), tolerance = 1e-12)

  quadrature <- .gk_single_panel(rule = 7L)
  first_node_times <- event_data$event_time[[1L]] * quadrature$nodes
  event_nodes <- event_data[rep(1L, length(first_node_times)), , drop = FALSE]
  event_nodes$event_time <- first_node_times
  expected_nodes <- .mm_event(standata$design_templates$event, event_nodes)
  expect_equal(unname(standata$W_gk[1L, , ]), unname(expected_nodes), tolerance = 1e-12)

  wrongly_scaled_nodes <- event_nodes
  wrongly_scaled_nodes$event_time <- first_node_times / standata$tmax
  expect_false(isTRUE(all.equal(
    unname(standata$W_gk[1L, , ]),
    unname(.mm_event(standata$design_templates$event, wrongly_scaled_nodes))
  )))

  expected_vcov_sd <- .mm(standata$design_templates$vcov$sd, event_data)
  expected_vcov_sd <- expected_vcov_sd[, colnames(expected_vcov_sd) != "(Intercept)", drop = FALSE]
  expect_equal(unname(standata$Xcov_sd), unname(expected_vcov_sd), tolerance = 1e-12)

  ordered_longitudinal <- longitudinal_data[
    order(longitudinal_data$id, longitudinal_data$marker, longitudinal_data$time),
    ,
    drop = FALSE
  ]
  expected_sigma <- .mm(
    standata$design_templates$distributional$sigma$default,
    ordered_longitudinal
  )
  expect_equal(dim(standata$X_sigma), dim(expected_sigma))
  expect_equal(as.numeric(standata$X_sigma), as.numeric(expected_sigma), tolerance = 1e-12)

  endpoint_baseline <- .mm(standata$basehaz_formula, event_data)
  expect_equal(
    unname(standata$Bs_event_c),
    unname(sweep(endpoint_baseline, 2L, standata$basehaz_col_means, "-")),
    tolerance = 1e-12
  )
})

test_that("time arrays retain their subject-by-ordinate order", {
  subject_data <- data.frame(id = c("a", "b"), time = 0, x = c(10, 20))
  times <- rbind(c(1, 2, 3), c(4, 5, 6))
  template <- .make_model_matrix_template(~ 0 + time + x, subject_data)

  evaluated <- .eval_template_on_times(template, subject_data, "time", times)

  time_column <- match("time", template$columns)
  x_column <- match("x", template$columns)
  expect_equal(evaluated[1L, , time_column], c(1, 2, 3))
  expect_equal(evaluated[2L, , time_column], c(4, 5, 6))
  expect_equal(evaluated[1L, , x_column], rep(10, 3))
  expect_equal(evaluated[2L, , x_column], rep(20, 3))
})

test_that("the Cox clock follows the Surv endpoint variable", {
  right_data <- data.frame(follow_up = c(1, 2), event = c(0, 1), x = 0)
  right <- .get_event_model_vars(
    survival::Surv(follow_up, event) ~ x,
    right_data
  )
  expect_identical(right$event_time_vars, "follow_up")
  expect_identical(right$event_stop_var, "follow_up")
  right_named <- .get_event_model_vars(
    survival::Surv(event = event, time = follow_up) ~ x,
    right_data
  )
  expect_identical(right_named$event_stop_var, "follow_up")

  counting_data <- data.frame(
    entry = c(0, 1),
    exit = c(1, 2),
    event = c(0, 1),
    x = 0
  )
  counting <- .get_event_model_vars(
    survival::Surv(entry, exit, event) ~ x,
    counting_data
  )
  expect_identical(counting$event_start_var, "entry")
  expect_identical(counting$event_stop_var, "exit")
  expect_identical(counting$event_time_vars, "exit")
})

test_that("simulate_joinme keeps spline-expanded random-effect column labels in truth draws", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time +
      (1 | id) +
      (1 + splines::bs(time, knots = 0.5, degree = 2) + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    n_id = 4,
    families = rep("gaussian", 2),
    time_cens = 1,
    times_obs = seq(0, 1, length.out = 3),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0))),
    seed = 31015,
    use_mirai = FALSE
  )

  expect_equal(
    colnames(sim$truth$re_draws$marker),
    c(
      "(Intercept)",
      "splines::bs(time, knots = 0.5, degree = 2)1",
      "splines::bs(time, knots = 0.5, degree = 2)2",
      "splines::bs(time, knots = 0.5, degree = 2)3"
    )
  )
  expect_equal(rownames(sim$truth$re_draws$marker), sim$marker_info$names)
})

test_that("fit and prediction support spline time terms in formulaLong", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = seq(0, 2, length.out = 4),
    quadrature_nodes = 7,
    seed = 3102,
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.15))),
  )

  formulaLong <- y ~ 1 + splines::ns(time, df = 3) + x1 +
    (0 + splines::bs(time, df = 3) | id) +
    (0 + x1 + (0 + splines::ns(time, df = 3) | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    formulaVCov = ~ 1,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      seed = 3102,
      refresh = 0
    )
  )

  expect_s3_class(fit, "JoiNMeFit")
  expect_false(any(grepl("^idx_time_", names(fit$stan_data))))

  pred <- posterior_epred(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    time_start = min(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 0.5, length.out = 15),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      n_samples = 10,
      chains = 1,
      iter_warmup = 5,
      iter_sampling = 5,
      refresh = 0
    ),
    seed = 3102
  )

  expect_s3_class(pred, "JoiNMeDynPred")
})
