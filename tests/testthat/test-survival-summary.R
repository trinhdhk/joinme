test_that("summary reports baseline survival covariates only when present", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(931)
  sim <- simulate_joinme(
    n_id = 6,
    families = rep("gaussian", 2),
    times_obs = seq(0, 3, length.out = 5),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2))),
    seed = 931
  )

  fit_cov <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 30,
      iter_sampling = 30,
      refresh = 0,
      seed = 931
    )
  )

  surv_tbl <- summary(fit_cov)$tables$survival_process
  expect_true(!is.null(surv_tbl))
  expect_true(all(c("x1", "x2") %in% surv_tbl$term))

  fit_none <- joinme(
    formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 30,
      iter_sampling = 30,
      refresh = 0,
      seed = 932
    )
  )

  expect_true(is.null(summary(fit_none)$tables$survival_process))
})
