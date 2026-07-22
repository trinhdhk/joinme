test_that("predict works with and without threading (cmdstanr)", {
  skip_on_cran()
#   skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  set.seed(101)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 101,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
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
    families = rep("gaussian", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 100,
      iter_sampling = 100,
      seed = 101,
      refresh = 0
    )
  )

  ndL <- sim$dataLong[sim$dataLong$id == 1, ]
  ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]

  pred1 <- suppressWarnings(
    posterior_epred(
      fit,
      newdataLong = ndL,
      newdataEvent = ndE,
      time_start = max(ndL$time),
      times = seq(max(ndL$time), max(ndL$time) + 0.5, length.out = 10),
      control = list(
        n_samples = 20,
        engine = "cmdstanr",
        chains = 1,
        iter_warmup = 50,
        iter_sampling = 1,
        threads_per_chain = 1,
        refresh = 0
      ),
      seed = 101
    )
  )
  expect_s3_class(pred1, "JoiNMeDynPred")

  pred2 <- suppressWarnings(
    posterior_epred(
      fit,
      newdataLong = ndL,
      newdataEvent = ndE,
      time_start = max(ndL$time),
      times = seq(max(ndL$time), max(ndL$time) + 0.5, length.out = 10),
      control = list(
        n_samples = 20,
        engine = "cmdstanr",
        chains = 1,
        iter_warmup = 50,
        iter_sampling = 1,
        threads_per_chain = 2,
        refresh = 0
      ),
      seed = 101
    )
  )
  expect_s3_class(pred2, "JoiNMeDynPred")
})
