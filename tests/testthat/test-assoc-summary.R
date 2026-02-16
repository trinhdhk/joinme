test_that("summary reports only active association components", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")

  set.seed(303)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 4, length.out = 8),
    seed = 303,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6)
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
    transforms = list(cv_total = list(type = "functional", expr = ~ softplus(x))),
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

  sum_obj <- summary(fit)
  assoc_tbl <- sum_obj$tables$assoc
  expect_true(any(assoc_tbl$term == "cv_total"))
  expect_false(any(assoc_tbl$term %in% c("cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker")))
  tf_tbl <- sum_obj$metadata$transform_formulas
  expect_true(any(tf_tbl$term == "cv_total"))
  expect_true(any(grepl("softplus", tf_tbl$formula)))
})
