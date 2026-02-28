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
  expect_true(any(assoc_tbl$term == "cv_total"))
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
