test_that("all", {
  skip_on_cran()
  set.seed(42)
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 6,
    D = 2,
    n_t = 4,
    seed = 42,
    include_marker_only = TRUE
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
    families = rep("student_t", 2),
    transforms = list(cv_total = list(type = "identity")),
    estimate_marker_weights = TRUE,
    marker_weight_scale = 0.5,
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 20,
      iter_sampling = 20,
      refresh = 0,
      adapt_delta = 0.95,
      max_treedepth = 12,
      seed = 42
    )
  )

  expect_s3_class(fit, "JoinMeFit")
  sum_obj <- summary(fit)
  expect_true(!is.null(sum_obj$tables))

  pred <- rstantools::posterior_epred(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    Tstart = max(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 1, length.out = 20),
    n_samples = 20,
    control = list(
      chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0
    )
  )

  expect_s3_class(pred, "JoinMeDynPred")
  p <- plot(pred, which = "survival", combined = FALSE)
  expect_true(inherits(p, "ggplot") || is.list(p))
})
