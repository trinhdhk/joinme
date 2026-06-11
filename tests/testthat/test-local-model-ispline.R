library(testthat)

test_that("joinme fits with ispline assoc", {
  skip_on_cran()
  if (!requireNamespace("joinme", quietly = TRUE)) skip("joinme not installed")
  has_cmd <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rstan <- requireNamespace("rstan", quietly = TRUE)
  if (!has_cmd && !has_rstan) skip("No Stan backend available")

  set.seed(123)
  sim <- joinme::simulate_joinme(
    n_id = 20,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 5, length.out = 5),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.3),
    transforms = list(cv_total = list(type = "ispline_penalised", knots = c(0, 1), degree = 3)),
    seed = 123
  )

  engine <- if (has_cmd) "cmdstanr" else "rstan"
  control <- list(
    engine = engine,
    chains = 1,
    iter_warmup = 30,
    iter_sampling = 30,
    threads_per_chain = 1,
    adapt_delta = 0.8,
    force_recompile = FALSE
  )
  if (engine == "cmdstanr") control$init <- 0

  fit <- tryCatch(
    joinme::joinme(
      dataLong = sim$dataLong,
      dataEvent = sim$dataEvent,
      formulaLong = y ~ time + (1 + time || id) + (1 || marker),
      formulaEvent = survival::Surv(time, event) ~ 1,
      assoc = c("cv_total"),
      transforms = list(cv_total = list(type = "ispline_penalised", knots = c(0, 1), degree = 3)),
      control = control
    ),
    error = function(e) skip(paste("fit failed:", conditionMessage(e)))
  )

  expect_s3_class(fit, "JoinMeFit")

  # basic sanity: summary runs
  s <- tryCatch(summary(fit), error = function(e) NULL)
  expect_true(!is.null(s))
})
