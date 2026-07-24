test_that("jm_family link specs are mapped into standata link_long", {
  set.seed(911)
  sim <- simulate_joinme(
    n_id = 8,
    families = c("gaussian", "bernoulli"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 3, length.out = 5),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2),
    seed = 911
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = list(
      jm_family("gaussian", link = "identity"),
      jm_family("bernoulli", link = "probit")
    )
  )

  expect_true(!is.null(sd$link_long))
  expect_equal(as.integer(sd$link_long), c(1L, 4L))
  expect_equal(as.integer(sd$inv_link_n_ops), c(1L, 2L))
  expect_equal(as.integer(sd$inv_link_ops[1, 1]), 0L)
  expect_equal(as.integer(sd$inv_link_ops[2, 1:2]), c(0L, 26L))

  sd2 <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = list(
      jm_family("gaussian", inv_link = ~ x),
      jm_family("bernoulli", inv_link = ~ inv_logit(x))
    )
  )
  expect_equal(as.integer(sd2$link_long), c(1L, 3L))
  expect_equal(as.integer(sd2$inv_link_ops[2, 1:2]), c(0L, 9L))

  sd3 <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = list(
      jm_family("gaussian", inv_link = ~ exp(x) + 0.1),
      jm_family("bernoulli", inv_link = ~ inv_logit(x))
    )
  )

  expect_equal(as.integer(sd3$link_long), c(0L, 3L))
  expect_true(is.na(sd3$link_names[1]))
  expect_equal(as.integer(sd3$inv_link_n_const[1]), 1L)
})

test_that("jm_family invertibility warning only triggers for likely non-injective inverse links", {
  expect_no_warning(
    jm_family("gaussian", inv_link = ~ exp(x) + 0.1)
  )

  expect_warning(
    jm_family("gaussian", inv_link = ~ abs(x)),
    "Custom inverse-link bytecode"
  )
})

test_that("extract.JoiNMeDynPred returns flattened draw payloads", {
  toy_draw <- matrix(rnorm(12), nrow = 3, ncol = 4)
  pred_obj <- JoiNMeDynPred$new(
    predictions = list(),
    quantiles = list(),
    draws = list(
      longitudinal = list(
        "1" = list(
          epred = list(matrix = toy_draw, marker_idx = c(1L, 1L, 2L, 2L), time = c(0, 1, 0, 1), scale = "epred")
        )
      ),
      survival = list(
        "1" = list(matrix = toy_draw, time = c(0, 1, 2, 3))
      )
    ),
    data = list(),
    metadata = list(scale = "epred", scales = c("epred")),
    call = quote(predict(fit)),
    tmax = 1,
    n_samples = 3
  )

  ext_long <- extract(pred_obj, what = "longitudinal", id = 1, scale = "epred")
  expect_true(is.list(ext_long$draws))
  expect_true("1" %in% names(ext_long$draws))
  expect_true("epred" %in% names(ext_long$draws[["1"]]))
  expect_equal(nrow(ext_long$draws[["1"]][["epred"]]), 3)

  ext_surv <- extract(pred_obj, what = "survival", id = 1)
  expect_true(is.list(ext_surv$draws))
  expect_equal(nrow(ext_surv$draws[["1"]]), 3)
})

test_that("extract.JoiNMeFit returns draw matrices by friendly names", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(707)
  sim <- simulate_joinme(
    n_id = 20,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.25),
    seed = 707
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
      iter_warmup = 60,
      iter_sampling = 60,
      refresh = 0,
      seed = 707
    )
  )

  ex_fixef <- extract(fit, what = "fixef")
  expect_true(is.matrix(ex_fixef$draws))
  expect_gt(nrow(ex_fixef$draws), 0)
  expect_gt(ncol(ex_fixef$draws), 0)

  ex_assoc <- extract(fit, what = "assoc")
  expect_true(is.matrix(ex_assoc$draws))
  expect_true(any(grepl("cv_total", colnames(ex_assoc$draws))))

  ex_eff <- extract(fit, what = "likelihood_scale", keep_chains = FALSE)
  expect_true(is.matrix(ex_eff$draws))
  expect_true(any(grepl("^beta_scaled: time$", colnames(ex_eff$draws))))
  expect_true(any(grepl("^id_sd_eff: time$", colnames(ex_eff$draws))))
  expect_true(any(grepl("^id_marker_row_scale_eff: time$", colnames(ex_eff$draws))))
})

test_that("extract.JoiNMeFit assoc uses canonical hazard-scale vcov coefficients", {
  draws_obj <- posterior::as_draws_matrix(stats::setNames(
    data.frame(
      raw = c(0.1, 0.1, 0.1),
      eff = c(0.4, 0.5, 0.6),
      check.names = FALSE
    ),
    c("z_alpha_vcov[1]", "alpha_vcov[1]")
  ))

  fit_obj <- structure(list(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(assoc_vcov = 1L, assoc_corr = 0L, assoc_cv_total = 0L, assoc_cv_mean = 0L,
                     assoc_cv_marker = 0L, assoc_cs_total = 0L, assoc_cs_mean = 0L, assoc_cs_marker = 0L,
                     D = 0L),
    config = list()
  ), class = "JoiNMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .package = "joinme"
  )

  ex_assoc <- extract(fit_obj, what = "assoc", keep_chains = FALSE)
  expect_equal(colnames(ex_assoc$draws), "vcov[1]")
  expect_equal(as.numeric(ex_assoc$draws[, 1]), c(0.4, 0.5, 0.6))
  expect_equal(as.character(ex_assoc$term_map$variable), "alpha_vcov[1]")
})

test_that("extract.JoiNMeFit fixef prefers internal coefficients on original-time scale", {
  draws_obj <- posterior::as_draws_matrix(stats::setNames(
    data.frame(
      beta_intercept = c(1.0, 1.0, 1.0),
      beta_time_raw = c(0.5, 0.5, 0.5),
      beta_scaled_intercept = c(1.0, 1.0, 1.0),
      beta_scaled_time = c(2.0, 2.0, 2.0),
      check.names = FALSE
    ),
    c("beta[1]", "beta[2]", "beta_scaled[1]", "beta_scaled[2]")
  ))

  fit_obj <- structure(list(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(P = 2L, x_cols = c("(Intercept)", "time")),
    config = list()
  ), class = "JoiNMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      if (is.null(variables)) return(draws_obj)
      posterior::subset_draws(draws_obj, variable = variables)
    },
    .package = "joinme"
  )

  ex_fixef <- extract(fit_obj, what = "fixef", keep_chains = FALSE)
  expect_equal(colnames(ex_fixef$draws), c("(Intercept)", "time"))
  expect_equal(as.numeric(ex_fixef$draws[, "time"]), c(0.5, 0.5, 0.5))
  expect_equal(
    as.character(ex_fixef$term_map$variable),
    c("beta[1]", "beta[2]")
  )
})

test_that("summary reports survival_process baseline covariates when present", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(1707)
  sim <- simulate_joinme(
    n_id = 20,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 6),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2),
    seed = 1707
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
      seed = 1707
    )
  )

  s <- summary(fit)
  expect_true(!is.null(s$tables$survival_process))
  surv_terms <- as.character(s$tables$survival_process$term)
  expect_true(all(c("x1", "x2") %in% surv_terms))
  expect_false(any(tolower(surv_terms) %in% c("(intercept)", "intercept", "1")))
})
