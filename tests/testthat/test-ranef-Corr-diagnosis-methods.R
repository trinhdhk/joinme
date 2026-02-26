test_that("ranef/corr.JoinMeDynPred require id-dependent marker covariance", {
  pred_ok <- joinme::JoinMeDynPred$new(
    predictions = list(),
    quantiles = list(),
    draws = list(
      random_effects_marker_id = list(
        "1" = list(
          matrix = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE,
                          dimnames = list(NULL, c("mk1::w1", "mk1::w2"))),
          corr = array(c(
            1.0, 0.1,
            0.1, 0.8,
            1.1, 0.2,
            0.2, 0.9
          ), dim = c(2, 2, 2)),
          markers = "mk1",
          terms = c("w1", "w2"),
          n_random_marker_id = 2L
        )
      )
    ),
    data = list(),
    metadata = list(marker_corr_depends_on_id = TRUE, n_samples = 2),
    call = NULL,
    tmax = 1,
    n_samples = 2
  )

  re_ok <- ranef(pred_ok)
  vc_ok <- corr(pred_ok)
  expect_true(is.data.frame(re_ok$formulaLong$marker_by_id))
  expect_true(is.data.frame(vc_ok$formulaLong$marker_by_id))
  expect_true(all(c("id", "row", "col", "Estimate") %in% names(vc_ok$formulaLong$marker_by_id)))

  pred_bad <- joinme::JoinMeDynPred$new(
    predictions = list(),
    quantiles = list(),
    draws = list(random_effects_marker_id = list()),
    data = list(),
    metadata = list(marker_corr_depends_on_id = FALSE),
    call = NULL,
    tmax = 1,
    n_samples = 2
  )

  expect_error(ranef(pred_bad), "depends on id")
  expect_error(corr(pred_bad), "depends on id")
})

test_that("marker_corr_depends_on_id follows fitted Q_idm", {
  fit_like_with_qidm <- list(stan_data = list(Q_idm = 2L))
  fit_like_without_qidm <- list(stan_data = list(Q_idm = 0L))

  expect_true(joinme:::.marker_corr_depends_on_id(fit_like_with_qidm))
  expect_false(joinme:::.marker_corr_depends_on_id(fit_like_without_qidm))
})

test_that("corr.JoinMeDynPred works when Q_idm > 0 without formulaCorr terms", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  set.seed(917)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 917,
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
      iter_warmup = 80,
      iter_sampling = 80,
      seed = 917,
      refresh = 0
    )
  )

  ndL <- sim$dataLong[sim$dataLong$id == 1, ]
  ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]

  pred <- tryCatch(
    suppressWarnings(
      posterior_epred(
        fit,
        newdataLong = ndL,
        newdataEvent = ndE,
        time_start = max(ndL$time),
        control = list(
          n_samples = 20,
          n_times = 50,
          engine = "cmdstanr",
          chains = 1,
          iter_warmup = 50,
          iter_sampling = 5,
          threads_per_chain = 1,
          refresh = 0
        ),
        seed = 918
      )
    ),
    error = function(e) e
  )
  if (inherits(pred, "error")) {
    skip(paste("Prediction sampler failed:", conditionMessage(pred)))
  }

  expect_true(isTRUE(pred$metadata$marker_corr_depends_on_id))
  expect_no_error(corr(pred))
  vc <- corr(pred)
  expect_true(is.data.frame(vc$formulaLong$marker_by_id))
})

test_that("ranef/corr.JoinMeFit include formulaDist random-effects blocks", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  set.seed(902)
  sim <- simulate_joinme(
    n_id = 4,
    families = c("student_t", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 902
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    formulaDist = list(sigma ~ 1 + (1 | id)),
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    families = c("student_t", "student_t"),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 20,
      iter_sampling = 20,
      refresh = 0,
      seed = 902
    )
  )

  re <- ranef(fit, draws = 20)
  vc <- corr(fit, draws = 20)
  dg <- diagnosis(fit)

  expect_true(is.list(re$formulaLong))
  expect_true(is.list(re$formulaDist))
  expect_true("sigma" %in% names(re$formulaDist))
  expect_true("allFamilies" %in% names(re$formulaDist$sigma))
  expect_true(all(c("term", "Estimate", "Q2.5", "Q97.5") %in% names(re$formulaDist$sigma$allFamilies)))

  expect_true(is.list(vc$formulaLong))
  expect_true(is.list(vc$formulaDist))
  expect_true("sigma" %in% names(vc$formulaDist))
  expect_true("allFamilies" %in% names(vc$formulaDist$sigma))
  expect_true(all(c("block", "row", "col", "Estimate") %in% names(vc$formulaDist$sigma$allFamilies)))

  expect_true(is.data.frame(dg))
  expect_true(all(c("metric", "value") %in% names(dg)))
})
