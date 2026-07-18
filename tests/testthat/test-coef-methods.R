make_mock_JoiNMe_fit_for_coef <- function(draws_obj, stan_data, config = list(), dataLong = NULL, dataEvent = NULL) {
  cache_env <- new.env(parent = emptyenv())
  fit_obj <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = stan_data,
      config = utils::modifyList(
        list(
          draws_default = NULL,
          transforms = list(),
          transforms_spec = list(),
          dist = list(dist_cols = list(), dist_re_terms = list())
        ),
        config
      ),
      dataLong = dataLong,
      dataEvent = dataEvent,
      cache_get = function(key) {
        if (!exists(key, envir = cache_env, inherits = FALSE)) return(NULL)
        get(key, envir = cache_env, inherits = FALSE)
      },
      cache_set = function(key, value) {
        assign(key, value, envir = cache_env)
        invisible(value)
      }
    ),
    class = "JoiNMeFit"
  )
  list(object = fit_obj, cache = cache_env)
}

test_that("posterior_fixef and posterior_ranef expose posterior extraction", {
  vars <- c(
    "beta[1]", "beta[2]",
    "u_id[1,1]", "u_id[1,2]", "u_id[2,1]", "u_id[2,2]",
    "v_marker[1,1]", "v_marker[1,2]",
    "w_idm[1,1,1]", "w_idm[2,1,1]",
    "gamma_w[1]",
    "beta_sigma[1]", "tau_sigma[1,1]", "z_sigma[1,1,1]",
    "alpha_L[1]", "lambda_L[1]", "z_L[1,1]", "z_L[2,1]"
  )
  vals <- c(
    rep(1.0, 4), rep(2.0, 4),
    rep(0.1, 4), rep(0.2, 4), rep(-0.1, 4), rep(-0.2, 4),
    rep(0.3, 4), rep(0.4, 4),
    rep(0.4, 4), rep(-0.4, 4),
    rep(0.5, 4),
    rep(0.7, 4), rep(0.5, 4), rep(0.2, 4),
    rep(0.6, 4), rep(0.5, 4), rep(0.4, 4), rep(-0.4, 4)
  )

  draws_obj <- posterior::as_draws_array(array(
    vals,
    dim = c(2, 2, length(vars)),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = vars
    )
  ))

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 2L,
      x_cols = c("(Intercept)", "time"),
      n_id = 2L,
      D = 1L,
      marker_levels = "m1",
      R_id = 2L,
      zid_cols = c("(Intercept)", "time"),
      R_mk = 2L,
      zmk_cols = c("(Intercept)", "time"),
      Q_idm = 1L,
      zidm_cols = "time",
      p_w = 1L,
      w_cols = "x1",
      K_event = 1L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      shared_marker_weights = 1L,
      n_re_sigma = 1L,
      K_sigma = 1L,
      G_sigma = 1L,
      K_cov = 0L,
      indep_idmarker_cov = 1L
    ),
    config = list(
      dist = list(
        dist_cols = list(sigma = "(Intercept)"),
        dist_re_terms = list(sigma = "(Intercept)")
      )
    ),
    dataLong = data.frame(id = c("id1", "id1", "id2", "id2")),
    dataEvent = data.frame(id = c("id1", "id2"))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = vars_keep)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      mat <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = vars_keep))
      as.matrix(mat[, vars_keep, drop = FALSE])
    },
    .package = "joinme"
  )

  fx <- fixef(fit_bundle$object, summary = FALSE)
  expect_true(is.matrix(fx))
  expect_equal(colnames(fx), c("(Intercept)", "time"))
  expect_equal(as.numeric(fx[, "time"]), rep(2, 4))
  expect_equal(posterior_fixef(fit_bundle$object), fx)

  re <- ranef(fit_bundle$object, summary = FALSE)
  expect_true(is.data.frame(re$formulaLong$id))
  expect_true(is.data.frame(re$formulaLong$marker))
  expect_true(is.data.frame(re$formulaLong$marker_by_id))
  expect_true(is.list(re$formulaDist$sigma))
  expect_true(is.data.frame(re$formulaDist$sigma$allFamilies))
  expect_true(is.data.frame(re$formulaVCov))
  expect_equal(unique(re$formulaLong$id$term), c("(Intercept)", "time"))
  expect_equal(unique(re$formulaLong$id$id), c("id1", "id2"))
  expect_equal(posterior_ranef(fit_bundle$object)$formulaLong$id$value, re$formulaLong$id$value)

  cf <- coef(fit_bundle$object, summary = FALSE)
  expect_true(is.data.frame(cf$formulaLong$population))
  expect_true(is.data.frame(cf$formulaLong$id))
  expect_true(is.data.frame(cf$formulaEvent))
  expect_true(is.list(cf$formulaDist$population))
  expect_true(is.list(cf$formulaDist$group_specific))
  expect_true(is.data.frame(cf$formulaVCov$id))

  id1_intercept <- cf$formulaLong$id$value[cf$formulaLong$id$id == "id1" & cf$formulaLong$id$term == "(Intercept)"]
  expect_equal(id1_intercept, rep(1.1, 4))
  id1_time <- cf$formulaLong$id$value[cf$formulaLong$id$id == "id1" & cf$formulaLong$id$term == "time"]
  expect_equal(id1_time, rep(2.2, 4))
  marker_time <- cf$formulaLong$marker$value[cf$formulaLong$marker$marker == "m1" & cf$formulaLong$marker$term == "time"]
  expect_equal(marker_time, rep(2.4, 4))
  marker_id_time <- cf$formulaLong$marker_by_id$value[cf$formulaLong$marker_by_id$id == "id1" & cf$formulaLong$marker_by_id$term == "time"]
  expect_equal(marker_id_time, rep(2.4, 4))

  sigma_combined <- cf$formulaDist$group_specific$sigma$allFamilies$value
  expect_equal(sigma_combined, rep(0.8, 4))
  vcov_id1 <- cf$formulaVCov$id$value[cf$formulaVCov$id$id == "id1"]
  expect_equal(vcov_id1, rep(0.2, 4))
})

test_that("coef summary reports grouped posterior summaries", {
  vars <- c("beta[1]", "u_id[1,1]")
  draws_obj <- posterior::as_draws_array(array(
    c(rep(1.0, 4), rep(0.2, 4)),
    dim = c(2, 2, 2),
    dimnames = list(iteration = c("1", "2"), chain = c("1", "2"), variable = vars)
  ))

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 1L,
      x_cols = "(Intercept)",
      n_id = 1L,
      D = 0L,
      R_id = 1L,
      zid_cols = "(Intercept)",
      R_mk = 0L,
      Q_idm = 0L,
      p_w = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L
    ),
    dataLong = data.frame(id = "id1")
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = vars_keep)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      mat <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = vars_keep))
      as.matrix(mat[, vars_keep, drop = FALSE])
    },
    .package = "joinme"
  )

  out <- coef(fit_bundle$object, summary = TRUE)
  expect_true(is.data.frame(out$formulaLong$id))
  expect_true(all(c("id", "term", "Estimate", "Q2.5", "Q97.5") %in% names(out$formulaLong$id)))
  expect_equal(out$formulaLong$id$Estimate, 1.2)
})

test_that("ranef and coef return original-time longitudinal terms", {
  vars <- c(
    "beta[1]", "beta[2]",
    "u_id[1,1]", "u_id[1,2]",
    "v_marker[1,1]", "v_marker[1,2]",
    "w_idm[1,1,1]", "w_idm[1,1,2]"
  )
  vals <- c(
    rep(1.0, 4), rep(0.5, 4),
    rep(0.2, 4), rep(4.0, 4),
    rep(0.3, 4), rep(6.0, 4),
    rep(0.4, 4), rep(8.0, 4)
  )

  draws_obj <- posterior::as_draws_array(array(
    vals,
    dim = c(2, 2, length(vars)),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = vars
    )
  ))

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 2L,
      x_cols = c("(Intercept)", "time"),
      n_id = 1L,
      D = 1L,
      marker_levels = "m1",
      R_id = 2L,
      zid_cols = c("(Intercept)", "time"),
      R_mk = 2L,
      zmk_cols = c("(Intercept)", "time"),
      Q_idm = 2L,
      zidm_cols = c("(Intercept)", "time"),
      p_w = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      tmax = 4,
      idx_time_uid = 2L,
      idx_time_vmk = 2L,
      idx_time_idm = 2L
    ),
    dataLong = data.frame(id = "id1"),
    dataEvent = data.frame(id = "id1")
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = vars_keep)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      mat <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = vars_keep))
      as.matrix(mat[, vars_keep, drop = FALSE])
    },
    .package = "joinme"
  )

  re_draws <- ranef(fit_bundle$object, summary = FALSE)
  expect_equal(unique(re_draws$formulaLong$id$value[re_draws$formulaLong$id$term == "time"]), 1)
  expect_equal(unique(re_draws$formulaLong$marker$value[re_draws$formulaLong$marker$term == "time"]), 1.5)
  expect_equal(unique(re_draws$formulaLong$marker_by_id$value[re_draws$formulaLong$marker_by_id$term == "time"]), 2)

  re_summary <- ranef(fit_bundle$object, summary = TRUE)
  expect_equal(re_summary$formulaLong$id$Estimate[re_summary$formulaLong$id$term == "time"], 1)
  expect_equal(re_summary$formulaLong$marker$Estimate[re_summary$formulaLong$marker$term == "time"], 1.5)

  cf <- coef(fit_bundle$object, summary = FALSE)
  expect_equal(unique(cf$formulaLong$id$value[cf$formulaLong$id$term == "time"]), 1.5)
  expect_equal(unique(cf$formulaLong$marker$value[cf$formulaLong$marker$term == "time"]), 2)
  expect_equal(unique(cf$formulaLong$marker_by_id$value[cf$formulaLong$marker_by_id$term == "time"]), 2.5)
})

test_that("coef extractors work on a fitted JoiNMe object", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(812)
  sim <- simulate_joinme(
    n_id = 16,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 5),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2),
    seed = 812
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 40,
      iter_sampling = 40,
      refresh = 0,
      seed = 812
    )
  )

  fx_summary <- fixef(fit)
  fx_draws <- posterior_fixef(fit)
  re_summary <- ranef(fit)
  re_draws <- posterior_ranef(fit)
  cf_summary <- coef(fit)
  cf_draws <- posterior_coef(fit)

  expect_true(is.data.frame(fx_summary))
  expect_true(is.matrix(fx_draws))
  expect_true(is.list(re_summary))
  expect_true(is.list(re_draws))
  expect_true(is.list(cf_summary))
  expect_true(is.list(cf_draws))

  expect_true(is.data.frame(cf_summary$formulaLong$population))
  expect_true(is.data.frame(cf_summary$formulaLong$id))
  expect_true(is.data.frame(cf_draws$formulaLong$population))
  expect_true(is.data.frame(cf_draws$formulaLong$id))
  expect_true(all(c("term", "Estimate", "Q2.5", "Q97.5") %in% names(cf_summary$formulaLong$population)))
  expect_true(all(c("draw", "id", "term", "value") %in% names(cf_draws$formulaLong$id)))
})

test_that("baseline hazard coefficients and predictions are extractable", {
  vars <- c("bs_gamma_c[1,1]", "bs_gamma_c[1,2]")
  vals <- c(
    0.1, 0.3,  # draw 1
    0.2, 0.4   # draw 2
  )
  draws_obj <- posterior::as_draws_array(array(
    vals,
    dim = c(2, 1, length(vars)),
    dimnames = list(iteration = c("1", "2"), chain = "1", variable = vars)
  ))

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 0L,
      p_w = 0L,
      K_event = 1L,
      Kbs = 2L,
      tmax = 2,
      Bs_event_c = matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE),
      basehaz = "formula",
      basehaz_cols = c("(Intercept)", "time")
    ),
    dataEvent = data.frame(id = c("id1", "id2"), time = c(1, 2))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = vars_keep)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      mat <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = vars_keep))
      as.matrix(mat[, vars_keep, drop = FALSE])
    },
    .get_draws_array = function(fit, variables = NULL, draws = NULL, seed = 1) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::as_draws_array(posterior::subset_draws(draws_obj, variable = vars_keep))
    },
    .package = "joinme"
  )

  bh_pred <- extract(fit_bundle$object, what = "basehaz", keep_chains = FALSE)
  expect_true(is.matrix(bh_pred$draws))
  expect_equal(nrow(bh_pred$draws), 2L)
  expect_equal(ncol(bh_pred$draws), 2L)
  expect_true(all(bh_pred$draws > 0))
  dm <- posterior::as_draws_matrix(draws_obj)
  coef_mat <- as.matrix(dm[, vars, drop = FALSE])
  expected <- exp(coef_mat %*% t(fit_bundle$object$stan_data$Bs_event_c)) / fit_bundle$object$stan_data$tmax
  expect_equal(unname(bh_pred$draws), unname(expected), tolerance = 1e-10)

  bh_draws <- draws(fit_bundle$object, what = "basehaz", format = "draws_matrix")
  expect_true(is.matrix(bh_draws))
  expect_true(all(grepl("id=", colnames(bh_draws), fixed = TRUE)))

  bh_pred_alias <- extract(fit_bundle$object, what = "baseline_hazard", keep_chains = FALSE)
  expect_equal(unname(bh_pred_alias$draws), unname(bh_pred$draws), tolerance = 1e-12)
})

test_that("summary reports intercept-only baseline hazard with intercept adjusted by log(tmax)", {
  vars <- c("bs_gamma_c[1,1]")
  vals <- c(log(0.2 * 2), log(0.3 * 2))
  draws_obj <- posterior::as_draws_array(array(
    vals,
    dim = c(2, 1, length(vars)),
    dimnames = list(iteration = c("1", "2"), chain = "1", variable = vars)
  ))

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 0L,
      p_w = 0L,
      family_long = 1L,
      family_names = "gaussian",
      n_id = 1L,
      R_id = 0L,
      R_mk = 0L,
      indep_id_re = 1L,
      indep_marker_re = 1L,
      K_event = 1L,
      Kbs = 1L,
      tmax = 2,
      Bs_event_c = matrix(1, nrow = 2, ncol = 1),
      basehaz = "formula",
      basehaz_cols = "(Intercept)"
    ),
    dataEvent = data.frame(id = c("id1", "id2"), time = c(1, 2))
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = vars_keep)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      mat <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = vars_keep))
      as.matrix(mat[, vars_keep, drop = FALSE])
    },
    .get_draws_array = function(fit, variables = NULL, draws = NULL, seed = 1) {
      if (is.null(variables)) return(draws_obj)
      vars_keep <- intersect(variables, posterior::variables(draws_obj))
      posterior::as_draws_array(posterior::subset_draws(draws_obj, variable = vars_keep))
    },
    .JoiNMe_sampler_diagnostics = function(fit) NULL,
    .package = "joinme"
  )

  sum_obj <- summary(fit_bundle$object, digits = 6)
  bh_tbl <- sum_obj$tables$baseline_hazard
  expect_true(is.data.frame(bh_tbl))
  expect_equal(bh_tbl$term, "(Intercept)")
  expect_true(all(c("Estimate", "Hazard.Ratio", "HR.Q2.5", "HR.Q97.5") %in% names(bh_tbl)))
  expected_log <- log(c(0.2, 0.3))
  expect_equal(bh_tbl$Estimate, mean(expected_log), tolerance = 1e-6)
  expect_equal(bh_tbl$Q2.5, as.numeric(stats::quantile(expected_log, 0.025)), tolerance = 1e-6)
  expect_equal(bh_tbl$Q97.5, as.numeric(stats::quantile(expected_log, 0.975)), tolerance = 1e-6)
  expect_equal(bh_tbl$Hazard.Ratio, exp(bh_tbl$Estimate), tolerance = 1e-6)
  expect_equal(bh_tbl$HR.Q2.5, exp(bh_tbl$Q2.5), tolerance = 1e-6)
  expect_equal(bh_tbl$HR.Q97.5, exp(bh_tbl$Q97.5), tolerance = 1e-6)
})

