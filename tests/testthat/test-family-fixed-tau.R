test_that("jm_family stores a valid fixed skew-Laplace quantile", {
  specification <- jm_family("skew_laplace", tau = 0.8)

  expect_s3_class(specification, "JoiNMe_family_spec")
  expect_identical(specification$family, "skew_double_exponential")
  expect_equal(specification$tau, 0.8)

  extracted <- .extract_family_and_link(specification)
  expect_identical(extracted$family_code, 9L)
  expect_equal(extracted$tau_fixed, 0.8)
})

test_that("jm_family enforces the quantile parameter's statistical support", {
  expect_error(
    jm_family("gaussian", tau = 0.5),
    "only defined for the skew-Laplace family"
  )
  expect_error(
    jm_family("skew_laplace", tau = 0),
    "strictly between zero and one"
  )
  expect_error(
    jm_family("skew_laplace", tau = 1),
    "not a valid skew-Laplace distribution"
  )
  expect_error(
    jm_family("skew_laplace", tau = c(0.2, 0.8)),
    "one finite numeric value"
  )
})

test_that("family validation preserves marker-specific fixed tau", {
  data_long <- data.frame(
    marker = factor(
      rep(c("gaussian_marker", "skew_fixed", "skew_estimated"), each = 2),
      levels = c("gaussian_marker", "skew_fixed", "skew_estimated")
    ),
    y = c(0.1, 0.2, -0.3, 0.4, 0.5, -0.2)
  )

  validated <- .validate_family_list(
    families = list(
      jm_family("gaussian"),
      jm_family("skew_laplace", tau = 0.8),
      jm_family("skew_laplace")
    ),
    D = 3L,
    dataLong = data_long,
    marker_var = "marker",
    y_var = "y"
  )

  expect_identical(validated$family_codes, c(1L, 9L, 9L))
  expect_identical(validated$use_tau_fixed, c(0L, 1L, 0L))
  expect_equal(validated$tau_fixed, c(0.5, 0.8, 0.5))

  tau_index <- .build_family_parameter_index(
    validated$family_codes,
    parameter = "tau",
    eligible = validated$use_tau_fixed == 0L
  )
  expect_identical(tau_index$n, 1L)
  expect_identical(tau_index$marker_to, c(0L, 0L, 1L))
})

test_that("standata carries marker-specific fixed tau into fitting", {
  simulation <- simulate_joinme(
    n_id = 5,
    families = c("gaussian", "skew_double_exponential"),
    times_obs = 0:2,
    seed = 991
  )

  standata <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = simulation$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = simulation$dataEvent,
    families = list(
      jm_family("gaussian"),
      jm_family("skew_laplace", tau = 0.8)
    ),
    assoc = "cv_total",
    seed = 991
  )

  expect_identical(standata$use_tau_fixed, c(0L, 1L))
  expect_equal(standata$tau_fixed, c(0.5, 0.8))
  expect_identical(standata$n_family_tau, 0L)
  expect_identical(standata$marker_to_tau_family, c(0L, 0L))
})

test_that("earlier scalar fixed-tau metadata remains usable in prediction", {
  normalised <- .normalise_fixed_tau_by_marker(
    use_tau_fixed = 1L,
    tau_fixed = 0.3,
    family_codes = c(1L, 9L, 9L)
  )

  expect_identical(normalised$use_tau_fixed, c(0L, 1L, 1L))
  expect_equal(normalised$tau_fixed, c(0.5, 0.3, 0.3))
})

test_that("distributional term mapping excludes fixed markers from estimated tau", {
  term_map <- .distributional_term_map(
    sd = list(
      D = 2L,
      family_long = c(9L, 9L),
      marker_levels = c("fixed_marker", "estimated_marker"),
      use_tau_fixed = c(1L, 0L),
      tau_fixed = c(0.8, 0.5),
      P_tau = 0L,
      n_re_tau = 0L
    ),
    cfg = list(),
    all_vars = "tau_family[1]"
  )

  expect_identical(
    unname(term_map[["tau_family[1]"]]),
    "tau_marker[estimated_marker]"
  )
})

test_that("summary reports fixed tau as a deterministic model quantity", {
  fit <- structure(
    list(
      stan_data = list(
        family_long = c(1L, 9L),
        marker_levels = c("gaussian_marker", "skew_marker"),
        use_tau_fixed = c(0L, 1L),
        tau_fixed = c(0.5, 0.8)
      ),
      config = list(family_long = c(1L, 9L))
    ),
    class = "JoiNMeFit"
  )

  summary_table <- .fixed_tau_summary_table(fit, digits = 3)
  expect_identical(summary_table$term, "tau_fixed[skew_marker]")
  expect_equal(summary_table$Estimate, 0.8)
  expect_equal(summary_table$Est.Error, 0)
  expect_equal(summary_table$Q2.5, 0.8)
  expect_equal(summary_table$Q97.5, 0.8)
  expect_true(is.na(summary_table$Rhat))
  expect_true(is.na(summary_table$ess_bulk))
  expect_true(is.na(summary_table$ess_tail))
})

test_that("summary.JoiNMeFit includes the family-fixed tau row", {
  raw_draws <- posterior::as_draws_array(array(
    c(
      0.1, 0.9, 1.1,
      0.2, 1.0, 1.2,
      0.0, 0.8, 1.0,
      0.3, 1.1, 1.3
    ),
    dim = c(2, 2, 3),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c("beta[1]", "sigma_family[1]", "sigma_family[2]")
    )
  ))
  fit <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 1L,
      p_w = 0L,
      K_event = 1L,
      Kbs = 0L,
      D = 2L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      x_cols = "(Intercept)",
      family_long = c(1L, 9L),
      family_names = c("gaussian", "skew_double_exponential"),
      marker_levels = c("gaussian_marker", "skew_marker"),
      use_tau_fixed = c(0L, 1L),
      tau_fixed = c(0.5, 0.8),
      P_sigma = 0L,
      P_nu = 0L,
      P_phi = 0L,
      P_alpha = 0L,
      P_kappa = 0L,
      P_tau = 0L,
      n_re_sigma = 0L,
      n_re_nu = 0L,
      n_re_phi = 0L,
      n_re_alpha = 0L,
      n_re_kappa = 0L,
      n_re_tau = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      marker_weight_sets_shared = 1L,
      tmax = 1
    ),
    formulaLong = y ~ 1 + (1 | id) + (1 | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(
      draws_default = NULL,
      family_long = c(1L, 9L),
      transform_spec = NULL,
      transforms = list(),
      dist = list(dist_cols = list())
    ),
    call = quote(joinme()),
    tmax = 1,
    dataLong = data.frame(
      y = c(0.1, 0.2),
      id = c(1, 2),
      marker = factor(c("gaussian_marker", "skew_marker"))
    ),
    dataEvent = data.frame(time = c(1, 1), event = c(1, 0))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      result <- raw_draws
      if (!is.null(variables)) {
        result <- posterior::subset_draws(result, variable = variables)
      }
      result
    },
    .joinme_sampler_diagnostics = function(fit) list(draws = 4),
    .package = "joinme"
  )

  fitted_summary <- summary(fit, include_corr = FALSE)
  fixed_row <- fitted_summary$tables$distributional[
    fitted_summary$tables$distributional$term == "tau_fixed[skew_marker]",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(fixed_row), 1L)
  expect_equal(fixed_row$Estimate, 0.8)
  expect_equal(fixed_row$Est.Error, 0)
})

test_that("fit, log-likelihood, and dynamic prediction use marker-indexed tau", {
  fit_data <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "longitudinal", "data", "fit.stan"),
    warn = FALSE
  )
  prediction_data <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "longitudinal", "data", "dynamic_prediction.stan"),
    warn = FALSE
  )
  fit_likelihood <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "etc", "functions", "joinme_fit_partial.stanfunctions"),
    warn = FALSE
  )
  prediction_likelihood <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "etc", "functions", "joinme_dynpred_partial.stanfunctions"),
    warn = FALSE
  )
  fit_log_likelihood <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "longitudinal", "generated_quantities", "fit_calculations.stan"),
    warn = FALSE
  )
  posterior_prediction <- readLines(
    testthat::test_path("..", "..", "inst", "stan", "include", "etc", "generated_quantities", "dynpred_output_calculations.stan"),
    warn = FALSE
  )

  expect_true(any(grepl("array\\[D\\].*use_tau_fixed", fit_data)))
  expect_true(any(grepl("tau_fixed[d]", fit_likelihood, fixed = TRUE)))
  expect_true(any(grepl("array\\[n_marker_types\\].*use_tau_fixed", prediction_data)))
  expect_true(any(grepl("tau_fixed[d]", prediction_likelihood, fixed = TRUE)))
  expect_true(any(grepl("tau_fixed[d]", fit_log_likelihood, fixed = TRUE)))
  expect_true(any(grepl("tau_fixed[d]", posterior_prediction, fixed = TRUE)))
})

test_that("simulation gives marker-specific family tau precedence", {
  fixed_simulation <- simulate_joinme(
    n_id = 4,
    families = jm_family("skew_laplace", tau = 0.8),
    family_params = list(
      skew_double_exponential = list(sigma = 1, tau = 0.3)
    ),
    times_obs = 0:1,
    seed = 733
  )

  expect_true(all(fixed_simulation$truth$distributional$rowwise$tau == 0.8))
  expect_identical(
    unname(fixed_simulation$truth$stan_fit$use_tau_fixed),
    1L
  )
  expect_equal(
    unname(fixed_simulation$truth$stan_fit$tau_fixed),
    0.8
  )
  expect_length(fixed_simulation$truth$stan_fit$tau_family, 0L)
  expect_identical(
    unname(fixed_simulation$truth$stan_fit$marker_to_tau_family),
    0L
  )

  shared_simulation <- simulate_joinme(
    n_id = 4,
    families = "skew_laplace",
    family_params = list(
      skew_double_exponential = list(sigma = 1, tau = 0.3)
    ),
    times_obs = 0:1,
    seed = 733
  )

  expect_true(all(shared_simulation$truth$distributional$rowwise$tau == 0.3))
  expect_identical(
    unname(shared_simulation$truth$stan_fit$use_tau_fixed),
    0L
  )
  expect_equal(
    unname(shared_simulation$truth$stan_fit$tau_family),
    0.3
  )
})

test_that("simulation supports fixed and modelled skew-Laplace tau by marker", {
  simulation <- simulate_joinme(
    n_id = 4,
    families = list(
      jm_family("skew_laplace", tau = 0.8),
      jm_family("skew_laplace")
    ),
    marker_levels = c("fixed_marker", "shared_marker"),
    family_params = list(
      skew_double_exponential = list(sigma = 1, tau = 0.3)
    ),
    formulaDist = list(tau ~ 1),
    truth = jm_truth(tau = list(
      intercept = c("(Intercept)" = stats::qlogis(0.3))
    )),
    times_obs = 0:1,
    seed = 734
  )

  rowwise <- simulation$truth$distributional$rowwise
  expect_true(all(rowwise$tau[rowwise$marker == "fixed_marker"] == 0.8))
  expect_equal(
    unique(rowwise$tau[rowwise$marker == "shared_marker"]),
    0.3
  )
  expect_identical(
    unname(simulation$truth$stan_fit$use_tau_fixed),
    c(1L, 0L)
  )
  expect_equal(
    unname(simulation$truth$stan_fit$tau_fixed),
    c(0.8, 0.5)
  )
  expect_equal(
    unname(simulation$truth$stan_fit$tau_family),
    0.3
  )
  expect_identical(
    unname(simulation$truth$stan_fit$marker_to_tau_family),
    c(0L, 1L)
  )
})

test_that("simulation rejects an unused tau regression", {
  expect_error(
    simulate_joinme(
      n_id = 4,
      families = jm_family("skew_laplace", tau = 0.8),
      formulaDist = list(tau ~ 1),
      times_obs = 0:1,
      seed = 735
    ),
    "has no simulated"
  )
})
