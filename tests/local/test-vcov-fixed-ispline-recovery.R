test_that("fixed vcov expit I-spline recovers a one-component association", {
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  expit_knots <- stats::plogis(c(-4, -2.5, -1, 0.5, 2, 4))
  tf <- joinme::joinme_tf(
    vcov = list(
      type = "ispline_expit",
      knots = expit_knots,
      coeff = c(0, 0.2, 0.45, 0.75, 0.95, 1.05, 1.1, 1.12),
      degree = 3
    )
  )

  sim <- joinme::simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time || id) + (0 + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    families = c("gaussian", "gaussian"),
    family_params = list(gaussian = list(sigma = 0.2)),
    n_id = 300,
    times_obs = seq(0, 10, length.out = 12),
    assoc = "vcov",
    truth = jm_truth(assoc_coef = list(slope = c(
      "vcov[1]" = 0.8, "vcov[2]" = 0, "vcov[3]" = 0
    ))),
    transforms = tf,
    seed = 2404,
    use_mirai = TRUE,
    n_workers = 2
  )

  fit <- joinme::joinme(
    formulaLong = y ~ 1 + time + (1 + time || id) + (0 + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    dataLong = sim$dataLong,
    dataEvent = sim$dataEvent,
    assoc = "vcov",
    transforms = tf,
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      iter_warmup = 250,
      iter_sampling = 250,
      adapt_delta = 0.95,
      max_treedepth = 12,
      refresh = 0,
      seed = 2404
    )
  )

  s <- suppressWarnings(summary(fit, draws = 200, seed = 2404))
  assoc_tbl <- subset(s$tables$assoc, grepl("^vcov\\[", term))

  expect_equal(nrow(assoc_tbl), 1L)
  expect_true(all(is.finite(assoc_tbl$Estimate)))
  expect_equal(assoc_tbl$Estimate[[1]], 0.8, tolerance = 0.35)
})
