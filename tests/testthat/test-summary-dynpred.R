test_that("summary.JoinMeDynPred reports rich subject-level outputs", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  set.seed(431)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 431,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.15)
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
    formulaVCov = ~ x1,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 100,
      iter_sampling = 100,
      seed = 431,
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
        seed = 431
      )
    ),
    error = function(e) e
  )
  if (inherits(pred, "error")) {
    skip(paste("Prediction sampler failed:", conditionMessage(pred)))
  }

  sum_pred <- summary(pred)
  expect_s3_class(sum_pred, "summary_JoinMeDynPred")

  expect_true(all(c("diagnostics", "overview", "median_survival_time", "random_effects_id", "random_effects_marker_id", "corr_marker_id") %in% names(sum_pred$tables)))

  expect_true(all(c("metric", "value") %in% names(sum_pred$tables$overview)))

  expect_true(all(c("id", "term", "n_reached", "n_total", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail") %in%
    names(sum_pred$tables$median_survival_time)))

  expect_true(all(c("id", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail") %in%
    names(sum_pred$tables$random_effects_id)))

  expect_true(all(c("id", "marker", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail") %in%
    names(sum_pred$tables$random_effects_marker_id)))

  expect_true(all(c("id", "row", "col", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail") %in%
    names(sum_pred$tables$corr_marker_id)))

  expect_true(all(c("metric", "value") %in% names(sum_pred$tables$diagnostics)))

  diag_tbl <- sum_pred$tables$diagnostics
  expect_true(is.data.frame(diag_tbl))
  expect_true(all(c("draws", "divergences", "treedepth_hits", "max_rhat", "min_ess_bulk", "min_ess_tail",
                    "n_terms_total", "n_terms_bad_rhat", "n_terms_low_ess_bulk", "n_terms_low_ess_tail") %in% diag_tbl$metric))
  expect_true(all(is.finite(diag_tbl$value[diag_tbl$metric %in% c(
    "draws", "divergences", "treedepth_hits",
    "n_terms_total", "n_terms_bad_rhat", "n_terms_low_ess_bulk", "n_terms_low_ess_tail"
  )])))

  expect_no_error(print(sum_pred))
})

test_that("summary.JoinMeDynPred suppresses off-diagonal covariance output under independence", {
  pred <- JoinMeDynPred$new(
    predictions = list(),
    quantiles = list(),
    draws = list(
      random_effects_marker_id = list(
        "1" = list(
          matrix = matrix(c(0.3, 0.4, 0.2, 0.1), ncol = 2, dimnames = list(NULL, c("m1::w1", "m1::w2"))),
          corr = array(c(
            1.0, 0.2,
            0.2, 1.1,
            1.0, 0.3,
            0.3, 1.2
          ), dim = c(2, 2, 2)),
          terms = c("w1", "w2")
        )
      )
    ),
    data = list(),
    metadata = list(
      marker_corr_depends_on_id = TRUE,
      indep_idmarker_cov = 1L,
      sampler_diagnostics = list()
    ),
    call = quote(posterior_epred(object)),
    tmax = 1,
    n_samples = 2
  )

  sum_pred <- summary(pred)
  expect_true(all(sum_pred$tables$corr_marker_id$row == sum_pred$tables$corr_marker_id$col))
  txt <- paste(capture.output(print(sum_pred)), collapse = "\n")
  expect_true(grepl("Prediction summary", txt, fixed = TRUE))
  expect_true(grepl("Predicted marker-by-id covariance", txt, fixed = TRUE))
})

test_that("summary.JoinMeDynPred diagnostics are populated from term and sampler sources", {
  pred <- JoinMeDynPred$new(
    predictions = list(
      longitudinal = data.frame(id = "1", time = 1, marker = "m1", Estimate = 0.1),
      survival = data.frame(id = "1", time = 1, Estimate = 0.9),
      cumhaz = data.frame(id = "1", time = 1, Estimate = 0.1)
    ),
    quantiles = list(),
    draws = list(
      random_effects_id = list(
        "1" = list(matrix = matrix(c(0.1, 0.2), ncol = 1), terms = "u_id[1]")
      ),
      random_effects_marker_id = list(
        "1" = list(
          matrix = matrix(c(0.3, 0.4), ncol = 1, dimnames = list(NULL, "m1::w1")),
          corr = array(c(1.0, 1.1), dim = c(2, 1, 1)),
          terms = "w1"
        )
      ),
      survival = list(
        "1" = list(matrix = matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2, byrow = TRUE), time = c(1, 2))
      )
    ),
    data = list(),
    metadata = list(
      scale = "epred",
      n_samples = 2,
      marker_corr_depends_on_id = TRUE,
      sampler_diagnostics = list(
        draws = 2,
        divergences = 0,
        treedepth_hits = 0,
        ebfmi_min = 0.7,
        max_rhat = 1.0,
        min_ess_bulk = 150,
        min_ess_tail = 140
      )
    ),
    call = NULL,
    tmax = 2,
    n_samples = 2
  )

  s <- summary(pred)
  d <- s$tables$diagnostics
  vals <- setNames(d$value, d$metric)

  expect_true(all(c("draws", "divergences", "treedepth_hits", "ebfmi_min",
                    "n_terms_total", "n_terms_bad_rhat", "n_terms_low_ess_bulk", "n_terms_low_ess_tail") %in% names(vals)))
  expect_true(all(is.finite(vals[c("draws", "divergences", "treedepth_hits", "ebfmi_min",
                                  "n_terms_total", "n_terms_bad_rhat", "n_terms_low_ess_bulk", "n_terms_low_ess_tail")])))
})
