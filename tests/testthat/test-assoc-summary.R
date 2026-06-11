test_that("summary reports only active association components", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")

  set.seed(303)
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 4,
    D = 2,
    n_t = 3,
    seed = 303,
    include_marker_only = TRUE
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
  expect_true(any(assoc_tbl$term == "cv_total (+)"))
  expect_false(any(assoc_tbl$term %in% c("cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker")))
})

test_that("summary hides marker weights when marker-weighted assoc terms are inactive", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  sim <- simulate_joinme(
    n_id = 40,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 5,
    times_obs = seq(0, 6, length.out = 9),
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0.2),
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
    fixed_marker_weights = FALSE,
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

test_that("summary prefers effective vcov association coefficients", {
  draws_obj <- posterior::as_draws_matrix(stats::setNames(
    data.frame(
      beta = c(0, 0, 0),
      raw = c(0.1, 0.1, 0.1),
      eff = c(0.4, 0.5, 0.6),
      check.names = FALSE
    ),
    c("beta[1]", "alpha_vcov[1]", "alpha_vcov_eff[1]")
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
  ), class = "JoinMeFit")

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

test_that("summary labels constrained weighted association terms and term-specific marker weights", {
  draws_obj <- posterior::as_draws_array(array(
    c(
      rep(0.5, 4),
      rep(0.4, 4),
      rep(1.0, 4),
      rep(-1.0, 4),
      rep(0.5, 4),
      rep(1.5, 4)
    ),
    dim = c(2, 2, 6),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c(
        "beta[1]",
        "alpha_cv_total",
        "alpha_cs_total",
        "marker_weights_eff_cv_total[1]",
        "marker_weights_eff_cs_total[1]",
        "marker_weights_eff_cs_total[2]"
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
      marker_levels = c("m1", "m2"),
      shared_marker_weights = 0L,
      Q_idm = 0L,
      indep_id_re = 0L,
      indep_marker_re = 0L,
      indep_idmarker_cov = 1L,
      family_names = "gaussian"
    ),
    config = list(draws_default = NULL),
    cache_get = function(key) NULL,
    cache_set = function(key, value) value
  ), class = "JoinMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .joinme_sampler_diagnostics = function(fit) list(),
    .package = "joinme"
  )

  assoc_tbl <- summary(fit_obj, include_corr = FALSE)$tables$assoc
  expect_true(any(assoc_tbl$term == "cv_total (+)"))
  expect_true(any(assoc_tbl$term == "cs_total (+)"))
  expect_true(any(grepl("^weight\\[cv_total\\]: m1$", assoc_tbl$term)))
  expect_true(any(grepl("^weight\\[cs_total\\]: m1$", assoc_tbl$term)))
  expect_true(any(grepl("^weight\\[cs_total\\]: m2$", assoc_tbl$term)))
})
