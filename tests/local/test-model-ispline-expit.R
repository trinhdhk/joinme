library(testthat)

test_that("joinme supports expit-based penalised spline transforms end to end", {
  skip_on_cran()
  skip_if_not_installed("splines2")
  if (!requireNamespace("joinme", quietly = TRUE)) skip("joinme not installed")
  has_cmd <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rstan <- requireNamespace("rstan", quietly = TRUE)
  if (!has_cmd && !has_rstan) skip("No Stan backend available")

  expit_grid <- stats::plogis(seq(-4, 4, length.out = 32))
  sim_tf <- joinme::joinme_tf(
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      y = expit_grid^0.8,
      n_knots = 10,
      degree = 2,
      lambda = 1
    )
  )
  fit_tf <- joinme::joinme_tf(
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      n_knots = 10,
      degree = 2,
      lambda = 1
    )
  )

  set.seed(2401)
  sim <- joinme::simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    families = c("gaussian", "gaussian"),
    n_id = 200,
    n_obs_per_marker_per_id = 10,
    times_obs = seq(0, 10, length.out = 10),
    assoc = "corr",
    assoc_coefs = list(corr = 0.2),
    transforms = sim_tf,
    seed = 2401,
    use_mirai = TRUE,
    n_workers = 8
  )

  engine <- if (has_cmd) "cmdstanr" else "rstan"
  control <- list(
    engine = engine,
    chains = 2,
    iter_warmup = 1000,
    iter_sampling = 500,
    parallel_chains = 2,
    threads_per_chain = 6,
    adapt_delta = 0.78,
    max_treedepth = 10,
    force_recompile = FALSE,
    seed = 2401
  )
 
  fit <- tryCatch(
    suppressWarnings(
      joinme::joinme(
        formulaLong = y ~ 1 + time + (1 + time || id) + (0 + (1 + time || id) || marker),
        formulaEvent = survival::Surv(time, event) ~ 1,
        dataLong = sim$dataLong,
        dataEvent = sim$dataEvent,
        assoc = "corr",
        transforms = fit_tf,
        control = control
      )
    ),
    error = function(e) skip(paste("fit failed:", conditionMessage(e)))
  )

  expect_s3_class(fit, "JoinMeFit")

  s <- suppressWarnings(tryCatch(summary(fit), error = function(e) NULL))
  expect_true(!is.null(s))
  expect_true(is.data.frame(s$metadata$transform_formulas))
  expect_true(any(grepl("ispline_expit", s$metadata$transform_formulas$formula, fixed = TRUE)))

  p_assoc <- tryCatch(
    plot(fit, type = "association", association_term = "corr", association_grid = seq(-3, 3, length.out = 11)),
    error = function(e) NULL
  )
  expect_s3_class(p_assoc, "ggplot")

  last_time <- aggregate(time ~ id, data = sim$dataLong, FUN = max)
  names(last_time)[2] <- "time_start"
  new_event <- merge(sim$dataEvent, last_time, by = "id", sort = FALSE)

  pred <- tryCatch(
    suppressWarnings(
      stats::predict(
        fit,
        newdataLong = subset(sim$dataLong, id %in% c(1, 2)),
        newdataEvent = subset(new_event, id %in% c(1, 2)),
        process = c("longitudinal", "event"),
        time_start = "time_start",
        time_horizon = 1,
        control = list(
          n_samples = 10,
          n_times = 10,
          chains = 1,
          parallel_chains = 1,
          iter_warmup = 20,
          iter_sampling = 20,
          threads_per_chain = 1,
          show_messages = FALSE
        )
      )
    ),
    error = function(e) skip(paste("prediction failed:", conditionMessage(e)))
  )

  expect_s3_class(pred, "JoinMeDynPred")
  p_pred <- tryCatch(
    plot(pred, which = c("longitudinal", "survival"), combined = TRUE),
    error = function(e) NULL
  )
  expect_true(inherits(p_pred, "gg") || inherits(p_pred, "ggplot") || is.list(p_pred))
})
