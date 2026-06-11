test_that("predict errors on marker levels not seen during fitting", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(551)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 4, length.out = 8),
    seed = 551,
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
    families = "student_t",
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 30,
      iter_sampling = 30,
      refresh = 0,
      seed = 551
    )
  )

  ndL <- sim$dataLong
  ndE <- sim$dataEvent

  # Inject an unseen marker level in prediction data
  ndL[["marker"]] <- as.character(ndL[["marker"]])
  ndL[["marker"]][1] <- "new_marker_level"
  ndL[["marker"]] <- factor(ndL[["marker"]])

  expect_error(
    posterior_epred(
      fit,
      newdataLong = ndL,
      newdataEvent = ndE,
      time_start = max(sim$dataLong$time),
      times = seq(0, max(sim$dataLong$time) + 0.5, length.out = 12),
      control = list(n_samples = 10),
      seed = 551
    ),
    regexp = "not seen during fitting|do not match fitted model levels|Cannot map marker value\\(s\\)"
  )
})
