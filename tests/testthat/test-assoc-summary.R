test_that("summary reports only active association components", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")

  set.seed(303)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("student_t", 2),
    times_obs = seq(0, 5, length.out = 3),
    censor_longitudinal_after_event = FALSE,
    seed = 303,
    use_mirai = FALSE
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("student_t", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0,
      seed = 303
    )
  )

  assoc_tbl <- summary(fit)$tables$assoc
  expect_true(any(assoc_tbl$term == "cv_total"))
  expect_false(any(assoc_tbl$term %in% c("cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker")))
})

test_that("summary hides marker weights when marker-weighted assoc terms are inactive", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  sim <- simulate_joinme(
    n_id = 40,
    families = rep("gaussian", 3),
    times_obs = seq(0, 6, length.out = 9),
    assoc = c("cv_mean"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_mean = 0.2))),
    seed = 908
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_mean"),
    families = rep("gaussian", 3),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = 2,
      iter_warmup = 80,
      iter_sampling = 80,
      refresh = 0,
      seed = 908
    )
  )

  assoc_tbl <- summary(fit)$tables$assoc
  expect_true(any(assoc_tbl$term == "cv_mean"))
  expect_false(any(grepl("^weight:", assoc_tbl$term)))
})

test_that("summary uses canonical hazard-scale vcov association coefficients", {
  draws_obj <- posterior::as_draws_matrix(stats::setNames(
    data.frame(
      beta = c(0, 0, 0),
      raw = c(0.1, 0.1, 0.1),
      eff = c(0.4, 0.5, 0.6),
      check.names = FALSE
    ),
    c("beta[1]", "z_alpha_vcov[1]", "alpha_vcov[1]")
  ))

  fit_obj <- structure(list(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 1L,
      x_cols = "(Intercept)",
      p_w = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 1L,
      D = 0L,
      Q_idm = 1L,
      indep_id_re = 0L,
      indep_marker_re = 0L,
      indep_idmarker_cov = 1L,
      family_names = "gaussian"
    ),
    config = list(draws_default = NULL),
    cache_get = function(key) NULL,
    cache_set = function(key, value) value
  ), class = "JoiNMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .joinme_sampler_diagnostics = function(fit) list(),
    .package = "joinme"
  )

  assoc_tbl <- summary(fit_obj, include_corr = FALSE)$tables$assoc
  expect_equal(assoc_tbl$term, "vcov[1]")
  expect_equal(assoc_tbl$Estimate, 0.5)
})

test_that("summary labels weighted association terms and term-specific marker weights", {
  draws_obj <- posterior::as_draws_array(array(
    c(
      rep(0.5, 4),
      rep(0.4, 4),
      rep(1.0, 4),
      rep(-1.0, 4),
      rep(0.0, 4),
      rep(0.5, 4),
      rep(1.5, 4),
      rep(0.25, 4),
      rep(-0.35, 4)
    ),
    dim = c(2, 2, 9),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c(
        "beta[1]",
        "alpha_cv_total",
        "alpha_cs_total",
        "marker_weights_eff_cv_total[1]",
        "marker_weights_eff_cv_total[2]",
        "marker_weights_eff_cs_total[1]",
        "marker_weights_eff_cs_total[2]",
        "marker_weight_mean[1]",
        "marker_weight_mean[2]"
      )
    )
  ))

  fit_obj <- structure(list(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 1L,
      x_cols = "(Intercept)",
      p_w = 0L,
      assoc_cv_total = 1L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 1L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      D = 2L,
      n_id = 0L,
      R_id = 0L,
      R_mk = 0L,
      marker_levels = c("m1", "m2"),
      marker_weight_offsets_by_term = list(
        cv_total = c(0, 0),
        cs_total = c(0, 0)
      ),
      marker_weight_sets_shared = 0L,
      n_marker_weight_means = 2L,
      marker_weight_set_cv_total = 1L,
      marker_weight_set_cs_total = 2L,
      Q_idm = 0L,
      indep_id_re = 0L,
      indep_marker_re = 0L,
      indep_idmarker_cov = 1L,
      family_names = "gaussian"
    ),
    config = list(
      draws_default = NULL,
      dist = list(dist_cols = list(), dist_re_terms = list())
    ),
    cache_get = function(key) NULL,
    cache_set = function(key, value) value
  ), class = "JoiNMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .joinme_sampler_diagnostics = function(fit) list(),
    .package = "joinme"
  )

  assoc_tbl <- summary(fit_obj, include_corr = FALSE)$tables$assoc
  expect_true(any(assoc_tbl$term == "cv_total"))
  expect_true(any(assoc_tbl$term == "cs_total"))
  expect_false(any(grepl("^weight", assoc_tbl$term)))

  weight_tbl <- summary(fit_obj, include_corr = FALSE)$tables$marker_weights
  expect_equal(
    weight_tbl$term,
    c("mean weight[cv_total]", "SD weight[cv_total]", "mean weight[cs_total]", "SD weight[cs_total]")
  )
  expect_equal(weight_tbl$Estimate[c(2L, 4L)], round(c(sqrt(0.8125), sqrt(2.0725)), 3))

  individual_weights <- marker_weights(fit_obj)
  expect_equal(individual_weights$set, rep(c("cv_total", "cs_total"), each = 2L))
  expect_equal(individual_weights$marker, rep(c("m1", "m2"), times = 2L))

  fixed_effect_weights <- fixef(fit_obj)
  expect_true(all(c("mean weight[cv_total]", "mean weight[cs_total]") %in% fixed_effect_weights$term))
  expect_false(any(grepl("^weight", fixed_effect_weights$term)))
  expect_equal(
    fixed_effect_weights$assoc_term[fixed_effect_weights$component == "marker_weight"],
    c("cv_total", "cs_total")
  )
  random_effect_weights <- ranef(fit_obj)
  expect_null(random_effect_weights$formulaLong$marker_weight)
  expect_equal(unique(random_effect_weights$assoc$assoc_term), c("cv_total", "cs_total"))
  expect_equal(random_effect_weights$assoc$marker, rep(c("m1", "m2"), 2L))
  combined_weights <- coef(fit_obj)
  expect_null(combined_weights$formulaLong$marker_weight)
  expect_equal(unique(combined_weights$assoc$assoc_term), c("cv_total", "cs_total"))
  expect_false("marker" %in% names(combined_weights$assoc))

  summary_object <- summary(fit_obj, include_corr = FALSE)
  expect_equal(summary_object$metadata$marker_weight_offsets$cv_total, c(m1 = 0, m2 = 0))
  expect_match(paste(capture.output(print(summary_object)), collapse = "\n"), "Marker-weight offsets \\[cv_total\\]")

  shared_object <- fit_obj
  shared_object$stan_data$marker_weight_sets_shared <- 1L
  shared_object$stan_data$n_marker_weight_means <- 1L
  shared_summary <- summary(shared_object, include_corr = FALSE)
  expect_equal(shared_summary$tables$marker_weights$term, c("mean weight", "SD weight"))
  shared_individual <- marker_weights(shared_object)
  expect_equal(shared_individual$set, rep("shared", 2L))
  expect_equal(shared_individual$marker, c("m1", "m2"))
})
