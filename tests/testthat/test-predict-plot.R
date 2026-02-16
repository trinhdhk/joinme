test_that("predict accepts time_start column and plot returns ribbons", {
  skip_on_cran()
  set.seed(2027)
  sim <- simulate_joinme(
    n_id = 6,
    families = rep("student_t", 3),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 6, length.out = 10),
    seed = 2027,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2),
    beta_event = c("(Intercept)" = -2.5, "x1" = 0.1, "x2" = -0.1),
    t_admin = 8
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time || id) +
    (0 + x1 + (1 + time || id) || marker)
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
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 100,
      iter_sampling = 50,
      refresh = 0,
      adapt_delta = 0.8,
      seed = 2027
    )
  )

  time_start <- tapply(sim$dataLong$time, sim$dataLong$id, max)
  ndE <- sim$dataEvent
  ndE$time_start <- as.numeric(time_start[as.character(ndE$id)])

  pred <- posterior_epred(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = ndE,
    time_start = "time_start",
    times = seq(0, max(sim$dataLong$time) + 1, length.out = 30),
    n_samples = 20,
    seed = 2027
  )

  expect_true(inherits(pred, "JoinMeDynPred"))
  expect_true(all(names(pred$metadata$conditioning_time_by_id) %in% as.character(ndE$id)))

  p_surv <- plot(pred, which = "survival", ci_type = "ribbon")
  expect_true(inherits(p_surv, "ggplot") || is.list(p_surv))
  if (inherits(p_surv, "ggplot")) {
    has_ribbon <- any(vapply(p_surv$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
    expect_true(has_ribbon)
  }

  if (requireNamespace("patchwork", quietly = TRUE) || requireNamespace("cowplot", quietly = TRUE)) {
    p_comb <- plot(pred, which = c("longitudinal", "survival", "cumhaz"), combined = TRUE)
    expect_true(inherits(p_comb, "gg") || is.list(p_comb))
  }

  tmp <- tempfile(fileext = ".png")
  p_to_save <- if (inherits(p_surv, "ggplot")) p_surv else p_surv[[1]]
  ggplot2::ggsave(tmp, p_to_save, width = 6, height = 4)
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})
