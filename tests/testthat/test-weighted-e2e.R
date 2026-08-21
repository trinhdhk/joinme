test_that("weighted grouping works end-to-end (simulate fit predict plot)", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  set.seed(2601)
  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | weighted(id, weights = id_w)) +
    (0 + x1 + (1 + time | weighted(id, weights = id_w)) | weighted(marker, weights = marker_w))
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sim <- simulate_joinme(
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    n_id = 5,
    families = rep("student_t", 2),
    times_obs = seq(0, 2, length.out = 6),
    seed = 2601,
    covariate_formulas = list(
      x1 ~ rnorm(n_id),
      x2 ~ rnorm(n_id),
      id_w ~ runif(n_id, 0.7, 1.4),
      marker_w ~ runif(n_id, 0.8, 1.6)
    ),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2))),
  )

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    formulaDist = list(
      sigma ~ 1 + (1 | weighted(id, weights = id_w))
    ),
    assoc = c("cv_total"),
    families = rep("student_t", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 60,
      iter_sampling = 60,
      refresh = 0,
      adapt_delta = 0.9,
      seed = 2601
    )
  )

  ndL <- sim$dataLong[sim$dataLong$id == sim$dataEvent$id[1], , drop = FALSE]
  ndE <- sim$dataEvent[sim$dataEvent$id == sim$dataEvent$id[1], , drop = FALSE]

  pred <- posterior_epred(
    fit,
    newdataLong = ndL,
    newdataEvent = ndE,
    time_start = max(ndL$time),
    seed = 2602,
    control = list(
      n_samples = 8,
      n_times = 50,
      engine = "cmdstanr",
      chains = 1,
      iter_warmup = 30,
      iter_sampling = 5,
      refresh = 0
    )
  )

  expect_s3_class(fit, "JoiNMeFit")
  expect_s3_class(pred, "JoiNMeDynPred")
  expect_true(all(is.finite(fit$stan_data$subject_weights)))
  expect_true(all(fit$stan_data$subject_weights > 0))

  p <- plot(pred, type = "survival", ci_type = "ribbon")
  expect_true(inherits(p, "ggplot") || is.list(p))
})
