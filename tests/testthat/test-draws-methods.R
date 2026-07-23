test_that("draws.JoiNMeFit returns cached renamed posterior draws", {
  raw_draws <- posterior::as_draws_array(array(
    c(
      0.2, 0.4, 0.6, -10,
      0.3, 0.5, 0.7, -9
    ),
    dim = c(2, 1, 4),
    dimnames = list(NULL, NULL, c("beta[1]", "beta[2]", "alpha_cv_total", "lp__"))
  ))

  fit_obj <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 2L,
      p_w = 0L,
      D = 0L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      x_cols = c("(Intercept)", "time"),
      assoc_cv_total = 1L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      shared_marker_weights = 1L
    ),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(draws_default = NULL, dist = list(dist_cols = list())),
    call = quote(joinme()),
    tmax = 1,
    dataLong = NULL,
    dataEvent = NULL
  )

  counter <- new.env(parent = emptyenv())
  counter$n <- 0L

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      counter$n <- counter$n + 1L
      out <- raw_draws
      if (!is.null(variables)) {
        out <- posterior::subset_draws(out, variable = variables)
      }
      out
    },
    .package = "joinme"
  )

  d1 <- draws(fit_obj, format = "draws_matrix")
  d2 <- draws(fit_obj, variables = "cv_", regex = TRUE, format = "draws_matrix")

  expect_equal(counter$n, 1L)
  expect_identical(colnames(d1), c("(Intercept)", "time", "cv_total", "lp__"))
  expect_identical(colnames(d2), "cv_total")

  arr <- as.array(fit_obj)
  expect_equal(dim(arr), c(2, 1, 4))
  expect_identical(dimnames(arr)[[3]], c("(Intercept)", "time", "cv_total", "lp__"))
})

test_that("summary and draws recover friendly terms from fits with stripped metadata", {
  raw_draws <- posterior::as_draws_array(array(
    c(
      0.2, 0.4, -0.3,
      0.3, 0.5, -0.2,
      0.1, 0.6, -0.1,
      0.4, 0.7, 0.0
    ),
    dim = c(2, 2, 3),
    dimnames = list(
      iteration = c("1", "2"),
      chain = c("1", "2"),
      variable = c("beta[1]", "beta[2]", "gamma_w[1]")
    )
  ))

  # This object reproduces fits created while the numerical Stan data were
  # stored without the character model-matrix metadata. The original formulas
  # and observed data remain sufficient to reconstruct the fitted term names.
  fit_obj <- JoiNMeFit$new(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 2L,
      p_w = 1L,
      K_event = 1L,
      Kbs = 0L,
      D = 0L,
      R_id = 0L,
      R_mk = 0L,
      Q_idm = 0L,
      assoc_cv_total = 0L,
      assoc_cv_mean = 0L,
      assoc_cv_marker = 0L,
      assoc_cs_total = 0L,
      assoc_cs_mean = 0L,
      assoc_cs_marker = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      shared_marker_weights = 1L,
      family_names = "gaussian",
      tmax = 3
    ),
    formulaLong = y ~ 1 + time + (1 | id) + (1 | marker),
    formulaEvent = survival::Surv(time, event) ~ x1,
    formulaVCov = NULL,
    config = list(
      draws_default = NULL,
      transform_spec = NULL,
      transforms = list(),
      dist = list(dist_cols = list())
    ),
    call = quote(joinme()),
    tmax = 3,
    dataLong = data.frame(
      y = c(0.1, 0.2, 0.3, 0.4),
      time = c(0, 1, 0, 1),
      id = c(1, 1, 2, 2),
      marker = factor(c("a", "b", "a", "b"))
    ),
    dataEvent = data.frame(
      time = c(2, 3),
      event = c(1, 0),
      x1 = c(-0.5, 0.5)
    )
  )

  # Emulate a stale cache produced by the former implementation. Versioned
  # cache names ensure that corrected labels are recomputed after package load.
  fit_obj$cache_set("renamed_draws_array", raw_draws)

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      out <- raw_draws
      if (!is.null(variables)) {
        out <- posterior::subset_draws(out, variable = variables)
      }
      out
    },
    .joinme_sampler_diagnostics = function(fit) list(draws = 4),
    .package = "joinme"
  )

  renamed <- draws(fit_obj, format = "draws_matrix")
  expect_identical(colnames(renamed), c("(Intercept)", "time", "x1"))

  fitted_summary <- summary(fit_obj, include_corr = FALSE)
  expect_identical(
    as.character(fitted_summary$tables$fixef$term),
    c("(Intercept)", "time")
  )
  expect_identical(
    as.character(fitted_summary$tables$survival_process$term),
    "x1"
  )
  expect_false(any(grepl("^(beta|gamma_w)\\[", unlist(lapply(
    fitted_summary$tables[c("fixef", "survival_process")],
    function(table) table$term
  )))))
})

test_that("draws.JoiNMeDynPred flattens stored prediction draws", {
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
      ),
      random_effects_id = list(
        "1" = list(matrix = matrix(c(0.1, 0.2, 0.3), ncol = 1), terms = "u_id[1]")
      )
    ),
    data = list(),
    metadata = list(scale = "epred", scales = c("epred")),
    call = quote(predict(fit)),
    tmax = 1,
    n_samples = 3
  )

  out <- draws(pred_obj, variables = c("survival", "random_effects_id"), regex = TRUE, format = "draws_matrix")
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 3)
  expect_true(any(grepl("survival", colnames(out), fixed = TRUE)))
  expect_true(any(grepl("random_effects_id", colnames(out), fixed = TRUE)))

  arr <- as.array(pred_obj)
  expect_equal(dim(arr)[1:2], c(3, 1))
  expect_gt(dim(arr)[3], 0)
})
