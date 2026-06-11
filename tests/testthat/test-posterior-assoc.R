make_mock_joinme_fit <- function(draws_obj, stan_data, config = list()) {
  cache_env <- new.env(parent = emptyenv())
  fit_obj <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = stan_data,
      config = utils::modifyList(
        list(draws_default = NULL, transforms = list(), transforms_spec = list(), dist = list(dist_cols = list())),
        config
      ),
      cache_get = function(key) {
        if (!exists(key, envir = cache_env, inherits = FALSE)) {
          return(NULL)
        }
        get(key, envir = cache_env, inherits = FALSE)
      },
      cache_set = function(key, value) {
        assign(key, value, envir = cache_env)
        invisible(value)
      }
    ),
    class = "JoinMeFit"
  )

  list(object = fit_obj, cache = cache_env)
}

mock_draws_array_subset <- function(draws_obj, variables = NULL) {
  arr <- if (is.null(variables)) {
    posterior::as_draws_array(draws_obj)
  } else {
    posterior::as_draws_array(posterior::subset_draws(draws_obj, variable = variables))
  }
  as.array(arr)
}

test_that("posterior_summary is an alias of summary for JoinMeFit", {
  draws_obj <- posterior::as_draws_array(array(
    c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6),
    dim = c(2, 2, 2),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c("beta[1]", "alpha_cv_mean")
    )
  ))

  fit_bundle <- make_mock_joinme_fit(
    draws_obj = draws_obj,
    stan_data = list(
      P = 1L,
      x_cols = "(Intercept)",
      p_w = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 1L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      D = 0L,
      Q_idm = 0L,
      R_id = 0L,
      R_mk = 0L,
      family_names = "gaussian"
    )
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .joinme_sampler_diagnostics = function(fit) list(),
    .package = "joinme"
  )

  summary_obj <- summary(fit_bundle$object, include_corr = FALSE)
  posterior_summary_obj <- posterior_summary(fit_bundle$object, include_corr = FALSE)

  expect_s3_class(summary_obj, "summary_JoinMeFit")
  expect_s3_class(posterior_summary_obj, "summary_JoinMeFit")
  expect_equal(posterior_summary_obj$tables$assoc, summary_obj$tables$assoc)
  expect_equal(posterior_summary_obj$tables$fixef, summary_obj$tables$fixef)
})

test_that("assoc propagates weighted association strength into marker-specific effects", {
  draws_obj <- posterior::as_draws_array(array(
    c(
      rep(2, 4),
      rep(1, 4),
      rep(3, 4)
    ),
    dim = c(2, 2, 3),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c("alpha_cv_total_eff", "marker_weights_eff[1]", "marker_weights_eff[2]")
    )
  ))

  fit_bundle <- make_mock_joinme_fit(
    draws_obj = draws_obj,
    stan_data = list(
      assoc_cv_total = 1L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      D = 2L,
      marker_levels = c("marker_a", "marker_b"),
      Q_idm = 0L,
      R_id = 0L,
      R_mk = 0L
    )
  )
  draw_array_calls <- 0L

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .get_draws_array = function(fit, variables, draws = NULL, seed = 1) {
      draw_array_calls <<- draw_array_calls + 1L
      mock_draws_array_subset(draws_obj, variables = variables)
    },
    .package = "joinme"
  )

  raw_assoc <- assoc(fit_bundle$object, summary = FALSE)
  expect_s3_class(raw_assoc, "PosteriorAssoc")
  expect_true(is.matrix(raw_assoc$cv_total))
  expect_equal(colnames(raw_assoc$cv_total), c("marker_a", "marker_b"))
  expect_true(all(raw_assoc$cv_total[, "marker_a"] == 1))
  expect_true(all(raw_assoc$cv_total[, "marker_b"] == 3))

  assoc_summary_obj <- assoc(fit_bundle$object, summary = TRUE)
  expect_s3_class(assoc_summary_obj, "PosteriorAssoc")
  expect_equal(assoc_summary_obj$cv_total$term, c("marker_a", "marker_b"))
  expect_equal(assoc_summary_obj$cv_total$Estimate, c(1, 3))
  expect_false(any(names(assoc_summary_obj) == "weight"))

  draw_array_calls_before_cache <- draw_array_calls
  cached_summary_obj <- assoc(fit_bundle$object, summary = TRUE)
  expect_equal(draw_array_calls, draw_array_calls_before_cache)
  expect_equal(cached_summary_obj$cv_total$Estimate, c(1, 3))

  alias_obj <- posterior_assoc(fit_bundle$object, summary = TRUE)
  expect_equal(alias_obj$cv_total$Estimate, assoc_summary_obj$cv_total$Estimate)

  printed <- paste(capture.output(print(assoc_summary_obj)), collapse = "\n")
  expect_true(grepl("Posterior association effects", printed, fixed = TRUE))
  expect_true(grepl("Association term: cv_total", printed, fixed = TRUE))
  expect_true(grepl("marker_a", printed, fixed = TRUE))
  expect_true(grepl("marker_b", printed, fixed = TRUE))
})

test_that("assoc uses different marker-weight structures for different weighted association terms", {
  draws_obj <- posterior::as_draws_array(array(
    c(
      rep(2, 4),
      rep(1, 4),
      rep(3, 4),
      rep(5, 4),
      rep(2, 4),
      rep(-2, 4)
    ),
    dim = c(2, 2, 6),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c(
        "alpha_cv_total_eff",
        "marker_weights_eff_cv_total[1]",
        "marker_weights_eff_cv_total[2]",
        "alpha_cs_total_eff",
        "marker_weights_eff_cs_total[1]",
        "marker_weights_eff_cs_total[2]"
      )
    )
  ))

  fit_bundle <- make_mock_joinme_fit(
    draws_obj = draws_obj,
    stan_data = list(
      assoc_cv_total = 1L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 1L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      D = 2L,
      marker_levels = c("marker_a", "marker_b"),
      shared_marker_weights = 0L,
      Q_idm = 0L,
      R_id = 0L,
      R_mk = 0L
    )
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .get_draws_array = function(fit, variables, draws = NULL, seed = 1) {
      mock_draws_array_subset(draws_obj, variables = variables)
    },
    .package = "joinme"
  )

  assoc_obj <- assoc(fit_bundle$object, summary = FALSE)
  expect_equal(as.numeric(assoc_obj$cv_total[, "marker_a"]), rep(1, 4))
  expect_equal(as.numeric(assoc_obj$cv_total[, "marker_b"]), rep(3, 4))
  expect_equal(as.numeric(assoc_obj$cs_total[, "marker_a"]), rep(5, 4))
  expect_equal(as.numeric(assoc_obj$cs_total[, "marker_b"]), rep(-5, 4))
})

test_that("summary labels covariance rows and columns with model terms", {
  draws_obj <- posterior::as_draws_array(array(
    c(0.1, 0.2, 0.3, 0.4),
    dim = c(2, 2, 1),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = "beta[1]"
    )
  ))

  fit_bundle <- make_mock_joinme_fit(
    draws_obj = draws_obj,
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
      assoc_vcov = 0L,
      D = 0L,
      Q_idm = 0L,
      R_id = 2L,
      R_mk = 1L,
      family_names = "gaussian",
      zid_cols = c("(Intercept)", "time"),
      zmk_cols = "marker_intercept"
    )
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .joinme_sampler_diagnostics = function(fit) list(),
    vcov.JoinMeFit = function(object, what = NULL, draws = NULL, ...) {
      if (identical(what, "id")) {
        return(data.frame(
          block = "id",
          row = c(1L, 2L, 2L),
          col = c(1L, 1L, 2L),
          Estimate = c(1.0, 0.2, 0.8),
          Est.Error = c(0.1, 0.1, 0.1),
          Q2.5 = c(0.8, 0.0, 0.6),
          Q97.5 = c(1.2, 0.4, 1.0),
          Rhat = c(1, 1, 1),
          ess_bulk = c(100, 100, 100),
          ess_tail = c(100, 100, 100),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        block = "marker",
        row = 1L,
        col = 1L,
        Estimate = 0.5,
        Est.Error = 0.1,
        Q2.5 = 0.3,
        Q97.5 = 0.7,
        Rhat = 1,
        ess_bulk = 100,
        ess_tail = 100,
        stringsAsFactors = FALSE
      )
    },
    .package = "joinme"
  )

  sum_obj <- summary(fit_bundle$object, include_corr = TRUE)

  expect_equal(sum_obj$tables$corr$id$row, c("(Intercept)", "time", "time"))
  expect_equal(sum_obj$tables$corr$id$col, c("(Intercept)", "(Intercept)", "time"))
  expect_equal(sum_obj$tables$corr$marker$row, "marker_intercept")
  expect_equal(sum_obj$tables$corr$marker$col, "marker_intercept")
})

test_that("assoc expands corr and vcov displays with explicit row and col labels", {
  draws_obj <- posterior::as_draws_array(array(
    c(
      rep(0.2, 4),
      rep(1.1, 4),
      rep(0.4, 4),
      rep(0.7, 4)
    ),
    dim = c(2, 2, 4),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c("alpha_corr_eff[1]", "alpha_vcov_eff[1]", "alpha_vcov_eff[2]", "alpha_vcov_eff[3]")
    )
  ))

  fit_bundle <- make_mock_joinme_fit(
    draws_obj = draws_obj,
    stan_data = list(
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 1L,
      assoc_vcov = 1L,
      D = 0L,
      Q_idm = 2L,
      R_id = 0L,
      R_mk = 0L,
      indep_idmarker_cov = 0L,
      zidm_cols = c("(Intercept)", "time")
    )
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) {
        return(draws_obj)
      }
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .get_draws_array = function(fit, variables, draws = NULL, seed = 1) {
      mock_draws_array_subset(draws_obj, variables = variables)
    },
    .package = "joinme"
  )

  assoc_summary_obj <- assoc(fit_bundle$object, summary = TRUE)
  expect_equal(names(assoc_summary_obj), c("corr", "vcov"))

  expect_equal(assoc_summary_obj$corr$component, 1L)
  expect_equal(assoc_summary_obj$corr$row, "time")
  expect_equal(assoc_summary_obj$corr$col, "(Intercept)")
  expect_equal(assoc_summary_obj$corr$term, "corr[time, (Intercept)]")
  expect_equal(assoc_summary_obj$corr$Estimate, 0.2)

  expect_equal(assoc_summary_obj$vcov$component, c(1L, 2L, 3L))
  expect_equal(assoc_summary_obj$vcov$row, c("(Intercept)", "time", "time"))
  expect_equal(assoc_summary_obj$vcov$col, c("(Intercept)", "(Intercept)", "time"))
  expect_equal(
    assoc_summary_obj$vcov$term,
    c("vcov[(Intercept), (Intercept)]", "vcov[time, (Intercept)]", "vcov[time, time]")
  )

  assoc_draws_obj <- assoc(fit_bundle$object, summary = FALSE)
  expect_equal(colnames(assoc_draws_obj$corr), "corr[time, (Intercept)]")
  expect_equal(
    colnames(assoc_draws_obj$vcov),
    c("vcov[(Intercept), (Intercept)]", "vcov[time, (Intercept)]", "vcov[time, time]")
  )

  printed <- paste(capture.output(print(assoc_summary_obj)), collapse = "\n")
  expect_true(grepl("Association term: corr", printed, fixed = TRUE))
  expect_true(grepl("Association term: vcov", printed, fixed = TRUE))
  expect_true(grepl("Displayed by matrix row using fitted random-effect term labels.", printed, fixed = TRUE))
  expect_true(grepl("row = time", printed, fixed = TRUE))
  expect_true(grepl("row = (Intercept)", printed, fixed = TRUE))
  expect_true(grepl("term", printed, fixed = TRUE))
  expect_true(grepl("Association term: corr", printed, fixed = TRUE))
  expect_true(grepl("Association term: vcov", printed, fixed = TRUE))
  expect_true(grepl("(Intercept)", printed, fixed = TRUE))
  expect_true(grepl("time", printed, fixed = TRUE))
})
