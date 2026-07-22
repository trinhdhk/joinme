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
