test_that("vcov association fits and predicts end-to-end", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    n_id = 4,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 3321,
    assoc = "vcov",
    assoc_coefs = list(vcov = c(0.2, 0.1, 0.25)),
    transforms = joinme_tf(vcov = "identity")
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    formulaVCov = ~ x1,
    assoc = "vcov",
    families = c("gaussian", "gaussian"),
    transforms = joinme_tf(vcov = "identity"),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 80,
      iter_sampling = 80,
      seed = 3321,
      refresh = 0
    )
  )

  fit_summary <- suppressWarnings(summary(fit, draws = 20, seed = 3321))
  expect_true(any(grepl("^vcov\\[", fit_summary$tables$assoc$term)))

  fit_assoc <- extract(fit, what = "assoc", keep_chains = FALSE)
  expect_true(any(grepl("^vcov\\[", colnames(fit_assoc$draws))))
  expect_true(any(grepl("^alpha_vcov_eff\\[", fit_assoc$term_map$variable)))
  expect_false(any(grepl("^alpha_vcov\\[", fit_assoc$term_map$variable)))

  ndL <- sim$dataLong[sim$dataLong$id == 1, , drop = FALSE]
  ndE <- sim$dataEvent[sim$dataEvent$id == 1, , drop = FALSE]
  pred <- suppressWarnings(
    posterior_epred(
      fit,
      newdataLong = ndL,
      newdataEvent = ndE,
      time_start = max(ndL$time),
      control = list(
        n_samples = 10,
        n_times = 25,
        engine = "cmdstanr",
        chains = 1,
        iter_warmup = 40,
        iter_sampling = 5,
        refresh = 0
      ),
      seed = 3321
    )
  )

  pred_summary <- summary(pred)
  expect_s3_class(pred_summary, "summary_JoinMeDynPred")
  expect_true("corr_marker_id" %in% names(pred_summary$tables))
})
