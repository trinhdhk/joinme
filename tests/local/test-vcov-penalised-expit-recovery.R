test_that("penalised expit vcov recovers a one-component association", {
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  expit_grid <- stats::plogis(seq(-4, 4, length.out = 25))
  sim_tf <- joinme::joinme_tf(
    vcov = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      y = expit_grid^0.8,
      n_knots = 8,
      degree = 2,
      lambda = 1
    )
  )
  fit_tf <- joinme::joinme_tf(
    vcov = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      n_knots = 8,
      degree = 2,
      lambda = 1
    )
  )

  sim <- joinme::simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time || id) + (0 + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    families = c("gaussian", "gaussian"),
    family_params = list(gaussian = list(sigma = 0.2)),
    n_id = 300,
    times_obs = seq(0, 10, length.out = 10),
    assoc = "vcov",
    truth = jm_truth(assoc_coef = list(slope = c(
      "vcov[1]" = 0.5, "vcov[2]" = 0, "vcov[3]" = 0
    ))),
    transforms = sim_tf,
    seed = 2405,
    use_mirai = TRUE,
    n_workers = 8
  )

  fit <- joinme::joinme(
    formulaLong = y ~ 1 + time + (1 + time || id) + (0 + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    dataLong = sim$dataLong,
    dataEvent = sim$dataEvent,
    assoc = "vcov",
    transforms = fit_tf,
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      iter_warmup = 1000,
      iter_sampling = 500,
      adapt_delta = 0.75,
      max_treedepth = 11,
      refresh = 200,
      threads_per_chain = 6,
      seed = 2405
    )
  )

  s <- suppressWarnings(summary(fit, draws = 200, seed = 2405))
  assoc_tbl <- subset(s$tables$assoc, grepl("^vcov\\[", term))

  expect_equal(nrow(assoc_tbl), 1L)
  expect_true(all(is.finite(assoc_tbl$Estimate)))
  expect_equal(assoc_tbl$Estimate[[1]], 0.5, tolerance = 0.2)
})
