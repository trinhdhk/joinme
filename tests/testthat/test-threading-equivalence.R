test_that("threaded and non-threaded Stan paths stay aligned", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Need >= 2 cores for threading test.")

  old_opts <- options(stan_preferred_engine = "cmdstanr")
  on.exit(options(old_opts), add = TRUE)

  seed <- 202406L
  sim <- simulate_joinme(
    n_id = 6,
    families = rep("student_t", 2),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 5, length.out = 10),
    seed = seed,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6),
    t_admin = 5.0,
    family_params = list(
      student_t = list(sigma = 0.35, nu = 4)
    )
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)

  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  ctrl_base <- list(
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 80,
    iter_sampling = 80,
    seed = seed,
    refresh = 0,
    adapt_delta = 0.8,
    max_treedepth = 12
  )

  fit_single <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    formulaCorr = ~1,
    assoc = "cv_total",
    transforms = list(
      cv_total = list(type = "identity"),
      corr = list(type = "identity")
    ),
    basehaz = "bs",
    n_knots = 3,
    basehaz_degree = 2,
    eps_fd = 1e-3,
    draws = 50,
    control = utils::modifyList(ctrl_base, list(threads_per_chain = 1L, grainsize = 1L))
  )

  fit_thread <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    formulaCorr = ~1,
    assoc = "cv_total",
    transforms = list(
      cv_total = list(type = "identity"),
      corr = list(type = "identity")
    ),
    basehaz = "bs",
    n_knots = 3,
    basehaz_degree = 2,
    eps_fd = 1e-3,
    draws = 50,
    control = utils::modifyList(ctrl_base, list(threads_per_chain = 2L, grainsize = 1L))
  )

  if (inherits(fit_thread$fit, "CmdStanMCMC")) {
    mod_thread <- tryCatch(fit_thread$fit$cmdstan_model(), error = function(e) NULL)
    threads_enabled <- tryCatch(joinme:::.cmdstan_threads_enabled(mod_thread), error = function(e) FALSE)
    if (!isTRUE(threads_enabled)) {
      skip("CmdStan model not compiled with threads; skipping threading equivalence checks.")
    }
  }

  vars <- c("beta[1]", "beta[2]", "alpha_cv_total")
  draws_single <- fit_single$fit$draws(variables = vars)
  draws_thread <- fit_thread$fit$draws(variables = vars)
  summ_single <- posterior::summarize_draws(draws_single, "mean", "sd")
  summ_thread <- posterior::summarize_draws(draws_thread, "mean", "sd")

  expect_true(all(vars %in% summ_single$variable))
  expect_true(all(vars %in% summ_thread$variable))

  means_single <- summ_single$mean[match(vars, summ_single$variable)]
  means_thread <- summ_thread$mean[match(vars, summ_thread$variable)]
  expect_equal(means_single, means_thread, tolerance = 0.2, scale = 1)

  pred_ctrl_base <- list(
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 1,
    seed = seed,
    refresh = 0
  )

  time_grid <- seq(0, max(sim$dataEvent$time) + 0.5, length.out = 15)

  pred_single <- predict(
    fit_single,
    sim$dataLong,
    sim$dataEvent,
    times = time_grid,
    control = utils::modifyList(pred_ctrl_base, list(n_samples = 30, n_times = length(time_grid), threads_per_chain = 1L, grainsize = 1L))
  )

  pred_thread <- predict(
    fit_thread,
    sim$dataLong,
    sim$dataEvent,
    times = time_grid,
    control = utils::modifyList(pred_ctrl_base, list(n_samples = 30, n_times = length(time_grid), threads_per_chain = 2L, grainsize = 1L))
  )

  surv_single <- pred_single$predictions$survival
  surv_thread <- pred_thread$predictions$survival
  expect_gt(nrow(surv_single), 0)
  expect_equal(nrow(surv_single), nrow(surv_thread))
  surv_merge <- merge(
    surv_single[, c("id", "time", "Survival")],
    surv_thread[, c("id", "time", "Survival")],
    by = c("id", "time"),
    suffixes = c("_single", "_thread")
  )
  expect_lt(max(abs(surv_merge$Survival_single - surv_merge$Survival_thread)), 0.1)

  long_single <- pred_single$predictions$longitudinal
  long_thread <- pred_thread$predictions$longitudinal
  expect_gt(nrow(long_single), 0)
  expect_equal(nrow(long_single), nrow(long_thread))
  long_merge <- merge(
    long_single[, c("id", "marker", "time", "Estimate")],
    long_thread[, c("id", "marker", "time", "Estimate")],
    by = c("id", "marker", "time"),
    suffixes = c("_single", "_thread")
  )
  expect_lt(max(abs(long_merge$Estimate_single - long_merge$Estimate_thread)), 0.15)

  # Check that plotting works (smoke test)
  expect_no_error(p_surv_s <- plot(pred_single, which = "survival"))
  expect_true(inherits(p_surv_s, "ggplot"))
  expect_no_error(p_long_s <- plot(pred_single, which = "longitudinal"))
  expect_true(inherits(p_long_s, "ggplot"))

  expect_no_error(p_surv_t <- plot(pred_thread, which = "survival"))
  expect_true(inherits(p_surv_t, "ggplot"))
  expect_no_error(p_long_t <- plot(pred_thread, which = "longitudinal"))
  expect_true(inherits(p_long_t, "ggplot"))
})
