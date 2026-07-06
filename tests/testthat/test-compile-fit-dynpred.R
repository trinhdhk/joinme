test_that("fit and dynpred compile with hardcoded GK nodes", {
  # Skip on CRAN to avoid long-running Stan compilation.
  skip_on_cran()

  # Require CmdStan for compilation checks.
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  # Simulate a small dataset to keep compilation and sampling lightweight.
  set.seed(2026)
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    quadrature_nodes = 7,
    seed = 2026,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  # Build longitudinal and survival formulas for a minimal joint model.
  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  # Fit with CmdStan to force compilation of the main model.
  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      seed = 2026,
      refresh = 0
    )
  )
  expect_s3_class(fit, "JoiNMeFit")

  # Run dynpred to compile the prediction program and exercise GK usage there.
  pred <- posterior_epred(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    time_start = min(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 0.5, length.out = 50),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      n_samples = 10,
      chains = 1,
      iter_warmup = 5,
      iter_sampling = 5,
      refresh = 0
    ),
    seed = 2026
  )
  expect_s3_class(pred, "JoiNMeDynPred")
})
