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
