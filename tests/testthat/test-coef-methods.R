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

test_that("coefficient methods expose summaries and draw-level extraction", {
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
      marker_weight_sets_shared = 1L,
      n_re_sigma = 1L,
      K_sigma = 1L,
      G_sigma = 1L,
      K_cov_sd = 0L,
      K_cov_corr = 0L,
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
  expect_equal(colnames(fx), c("(Intercept)", "time", "event: x1"))
  expect_equal(as.numeric(fx[, "time"]), rep(2, 4))
  expect_equal(fixef(fit_bundle$object, summary = FALSE), fx)
  expect_equal(extract(fit_bundle$object, what = "fixed_effects")$posterior_draws, fx)
  expect_equal(posterior_summary(fit_bundle$object, what = "fixef", summary = FALSE), fx)

  re <- ranef(fit_bundle$object, summary = FALSE)
  expect_true(is.data.frame(re$formulaLong$id))
  expect_true(is.data.frame(re$formulaLong$marker))
  expect_true(is.data.frame(re$formulaLong$marker_by_id))
  expect_true(is.list(re$formulaDist$sigma))
  expect_true(is.data.frame(re$formulaDist$sigma$allFamilies))
  expect_true(is.data.frame(re$formulaVCov))
  expect_equal(unique(re$formulaLong$id$term), c("(Intercept)", "time"))
  expect_equal(unique(re$formulaLong$id$id), c("id1", "id2"))
  expect_equal(ranef(fit_bundle$object, summary = FALSE)$formulaLong$id$value, re$formulaLong$id$value)
  expect_equal(extract(fit_bundle$object, what = "random_effects")$posterior_draws, re)
  expect_equal(posterior_summary(fit_bundle$object, what = "ranef", summary = FALSE), re)

  cf <- coef(fit_bundle$object, summary = FALSE)
  expect_null(cf$formulaLong$population)
  expect_true(is.data.frame(cf$formulaLong$id))
  expect_true(is.data.frame(cf$formulaEvent))
  expect_true(is.list(cf$formulaDist$population))
  expect_true(is.list(cf$formulaDist$group_specific))
  expect_true(is.data.frame(cf$formulaVCov$id))
  expect_equal(extract(fit_bundle$object, what = "coefficients")$posterior_draws, cf)
  expect_equal(posterior_summary(fit_bundle$object, what = "coef", summary = FALSE), cf)

  fixed_intervals <- posterior_interval(fit_bundle$object, prob = 0.8, what = "fixef")
  expect_true(is.matrix(fixed_intervals))
  expect_identical(colnames(fixed_intervals), c("10%", "90%"))
  expect_identical(rownames(fixed_intervals), c("(Intercept)", "time", "event: x1"))
  expect_equal(unname(fixed_intervals["time", ]), c(2, 2))

  random_intervals <- posterior_interval(fit_bundle$object, prob = 0.8, what = "ranef")
  expect_true(is.data.frame(random_intervals$formulaLong$id))
  expect_true(all(c("10%", "90%") %in% names(random_intervals$formulaLong$id)))
  expect_false(any(c("draw", "value", "fixed", "random") %in% names(random_intervals$formulaLong$id)))

  coefficient_intervals <- posterior_interval(fit_bundle$object, prob = 0.8, what = "coef")
  expect_true(is.data.frame(coefficient_intervals$formulaLong$id))
  expect_true(all(c("10%", "90%") %in% names(coefficient_intervals$formulaLong$id)))
  expect_equal(
    coefficient_intervals$formulaLong$id$`10%`[
      coefficient_intervals$formulaLong$id$id == "id1" &
        coefficient_intervals$formulaLong$id$term == "(Intercept)"
    ],
    1.1
  )

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

test_that("removed posterior coefficient entry points are absent", {
  package_namespace <- asNamespace("joinme")
  removed_entry_points <- paste0("posterior_", c("fixef", "ranef", "coef"))
  expect_false(any(vapply(
    removed_entry_points,
    exists,
    logical(1),
    envir = package_namespace,
    inherits = FALSE
  )))
})

test_that("marker weights follow fixed, random, and combined coefficient semantics", {
  variables <- c(
    "beta[1]",
    "gamma_w[1]",
    "marker_weight_mean[1]", "marker_weight_mean[2]",
    "marker_weights_eff_cv_total[1]", "marker_weights_eff_cv_total[2]",
    "marker_weights_eff_cs_marker[1]", "marker_weights_eff_cs_marker[2]"
  ) # population coefficient, two weight-set means, and four effective marker weights
  values <- c(
    rep(1.0, 4),
    rep(0.7, 4),
    rep(-0.2, 4), rep(0.5, 4),
    rep(1.2, 4), rep(0.1, 4),
    rep(0.0, 4), rep(1.1, 4)
  ) # constant draws chosen so the centred departures have simple exact values
  draws_obj <- posterior::as_draws_array(array(
    values,
    dim = c(2, 2, length(variables)),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = variables
    )
  )) # two-chain posterior fixture retaining valid array dimensions

  fit_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = list(
      P = 1L,
      x_cols = "(Intercept)",
      n_id = 0L,
      D = 2L,
      marker_levels = c("m1", "m2"),
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      p_w = 1L,
      w_cols = "treatment",
      K_event = 1L,
      assoc_cv_total = 1L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 1L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      marker_weight_sets_shared = 0L,
      n_marker_weight_means = 2L,
      marker_weight_set_cv_total = 2L,
      marker_weight_set_cs_marker = 1L,
      marker_weight_offsets_by_term = list(
        cv_total = c(0.2, -0.4),
        cs_marker = c(0.3, 0.8)
      )
    )
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      variables_kept <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = variables_kept)
    },
    .get_draws_matrix = function(fit, variables = NULL, draws = NULL, seed = 1) {
      variables_kept <- intersect(variables, posterior::variables(draws_obj))
      matrix_draws <- posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = variables_kept))
      as.matrix(matrix_draws[, variables_kept, drop = FALSE])
    },
    .get_draws_array = function(fit, variables = NULL, draws = NULL, seed = 1) {
      variables_kept <- intersect(variables, posterior::variables(draws_obj))
      posterior::subset_draws(draws_obj, variable = variables_kept)
    },
    .package = "joinme"
  )

  fixed_draws <- fixef(fit_bundle$object, summary = FALSE)
  expect_equal(unname(colnames(fixed_draws)), c("(Intercept)", "event: treatment", "mean weight[cs_marker]", "mean weight[cv_total]"))
  expect_false(any(grepl("^weight", colnames(fixed_draws))))

  fixed_summary <- fixef(fit_bundle$object, summary = TRUE)
  expect_true(all(c("component", "event", "assoc_term", "term") %in% names(fixed_summary)))
  expect_equal(fixed_summary$component, c("longitudinal", "event", "marker_weight", "marker_weight"))
  expect_equal(fixed_summary$event[fixed_summary$component == "event"], "event")
  expect_equal(fixed_summary$term[fixed_summary$component == "event"], "treatment")
  expect_equal(
    fixed_summary$assoc_term[fixed_summary$component == "marker_weight"],
    c("cs_marker", "cv_total")
  )

  random_draws <- ranef(fit_bundle$object, summary = FALSE)
  expect_null(random_draws$formulaLong$marker_weight)
  expect_true(is.data.frame(random_draws$assoc))
  expect_true(all(c("assoc_term", "marker", "term", "value") %in% names(random_draws$assoc)))
  departure_values <- random_draws$assoc
  expect_equal(unique(departure_values$value[departure_values$assoc_term == "cv_total" & departure_values$marker == "m1"]), 0.5)
  expect_equal(unique(departure_values$value[departure_values$assoc_term == "cv_total" & departure_values$marker == "m2"]), 0)
  expect_equal(unique(departure_values$value[departure_values$assoc_term == "cs_marker" & departure_values$marker == "m1"]), -0.1)
  expect_equal(unique(departure_values$value[departure_values$assoc_term == "cs_marker" & departure_values$marker == "m2"]), 0.5)

  combined_draws <- coef(fit_bundle$object, summary = FALSE)
  expect_null(combined_draws$formulaLong$marker_weight)
  expect_true(is.data.frame(combined_draws$assoc))
  expect_false("marker" %in% names(combined_draws$assoc))
  effective_values <- combined_draws$assoc
  expect_equal(unique(effective_values$value[effective_values$assoc_term == "cv_total" & effective_values$term == "m1"]), 1.2)
  expect_equal(unique(effective_values$value[effective_values$assoc_term == "cs_marker" & effective_values$term == "m2"]), 1.1)
  expect_null(combined_draws$formulaLong$population)
  expect_equal(unique(combined_draws$formulaEvent$term), "treatment")

  random_summary <- ranef(fit_bundle$object, summary = TRUE)$assoc
  expect_true(all(c("assoc_term", "marker", "Estimate") %in% names(random_summary)))
  expect_equal(random_summary$Estimate[random_summary$assoc_term == "cv_total" & random_summary$marker == "m1"], 0.5)

  combined_summary <- coef(fit_bundle$object, summary = TRUE)$assoc
  expect_true(all(c("assoc_term", "term", "Estimate") %in% names(combined_summary)))
  expect_false("marker" %in% names(combined_summary))
  expect_equal(combined_summary$Estimate[combined_summary$assoc_term == "cs_marker" & combined_summary$term == "m2"], 1.1)

  shared_object <- fit_bundle$object
  shared_object$stan_data$marker_weight_sets_shared <- 1L
  shared_object$stan_data$n_marker_weight_means <- 1L
  shared_object$stan_data$marker_weight_set_cs_marker <- 1L
  shared_weights <- .marker_weight_long_draws(shared_object, quantity = "effective")
  expect_equal(unique(shared_weights$assoc_term), "shared")
  expect_equal(unique(shared_weights$marker), c("m1", "m2"))
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
      tmax = 4
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
  expect_equal(unique(re_draws$formulaLong$id$value[re_draws$formulaLong$id$term == "time"]), 4)
  expect_equal(unique(re_draws$formulaLong$marker$value[re_draws$formulaLong$marker$term == "time"]), 6)
  expect_equal(unique(re_draws$formulaLong$marker_by_id$value[re_draws$formulaLong$marker_by_id$term == "time"]), 8)

  re_summary <- ranef(fit_bundle$object, summary = TRUE)
  expect_equal(re_summary$formulaLong$id$Estimate[re_summary$formulaLong$id$term == "time"], 4)
  expect_equal(re_summary$formulaLong$marker$Estimate[re_summary$formulaLong$marker$term == "time"], 6)

  cf <- coef(fit_bundle$object, summary = FALSE)
  expect_equal(unique(cf$formulaLong$id$value[cf$formulaLong$id$term == "time"]), 4.5)
  expect_equal(unique(cf$formulaLong$marker$value[cf$formulaLong$marker$term == "time"]), 6.5)
  expect_equal(unique(cf$formulaLong$marker_by_id$value[cf$formulaLong$marker_by_id$term == "time"]), 8.5)
})

test_that("coef extractors work on a fitted JoiNMe object", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(812)
  sim <- simulate_joinme(
    n_id = 16,
    families = rep("gaussian", 2),
    times_obs = seq(0, 4, length.out = 5),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2))),
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
  fx_draws <- fixef(fit, summary = FALSE)
  re_summary <- ranef(fit)
  re_draws <- ranef(fit, summary = FALSE)
  cf_summary <- coef(fit)
  cf_draws <- coef(fit, summary = FALSE)

  expect_true(is.data.frame(fx_summary))
  expect_true(is.matrix(fx_draws))
  expect_true(is.list(re_summary))
  expect_true(is.list(re_draws))
  expect_true(is.list(cf_summary))
  expect_true(is.list(cf_draws))

  expect_null(cf_summary$formulaLong$population)
  expect_true(is.data.frame(cf_summary$formulaLong$id))
  expect_null(cf_draws$formulaLong$population)
  expect_true(is.data.frame(cf_draws$formulaLong$id))
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
  expect_true(is.matrix(bh_pred$posterior_draws))
  expect_equal(nrow(bh_pred$posterior_draws), 2L)
  expect_equal(ncol(bh_pred$posterior_draws), 2L)
  expect_true(all(bh_pred$posterior_draws > 0))
  dm <- posterior::as_draws_matrix(draws_obj)
  coef_mat <- as.matrix(dm[, vars, drop = FALSE])
  expected <- exp(coef_mat %*% t(fit_bundle$object$stan_data$Bs_event_c)) / fit_bundle$object$stan_data$tmax
  expect_equal(unname(bh_pred$posterior_draws), unname(expected), tolerance = 1e-10)

  bh_draws <- posterior_draws(fit_bundle$object, what = "basehaz", format = "draws_matrix")
  expect_true(is.matrix(bh_draws))
  expect_true(all(grepl("id=", colnames(bh_draws), fixed = TRUE)))

  bh_pred_alias <- extract(fit_bundle$object, what = "baseline_hazard", keep_chains = FALSE)
  expect_equal(unname(bh_pred_alias$posterior_draws), unname(bh_pred$posterior_draws), tolerance = 1e-12)
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
    .joinme_sampler_diagnostics = function(fit) NULL,
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

  spline_bundle <- make_mock_JoiNMe_fit_for_coef(
    draws_obj = draws_obj,
    stan_data = utils::modifyList(
      fit_bundle$object$stan_data,
      list(basehaz = "bs", basehaz_cols = "1")
    ),
    dataEvent = data.frame(id = c("id1", "id2"), time = c(1, 2))
  )
  spline_summary <- summary(spline_bundle$object, digits = 6)
  expect_equal(spline_summary$tables$baseline_hazard$term, "basis_1")
})

