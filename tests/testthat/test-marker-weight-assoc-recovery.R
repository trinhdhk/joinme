test_that("fit recovers signed marker weights and alpha without passing marker_weights", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    n_id = 120,
    families = rep("gaussian", 3),
    n_obs_per_marker_per_id = 8,
    times_obs = seq(0, 8, length.out = 12),
    marker_weights = c(0.5, 1, -0.5),
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.6),
    transforms = list(cv_total = list(type = "identity")),
    seed = 4102
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = rep("gaussian", 3),
    transforms = list(cv_total = list(type = "identity")),
    estimate_marker_weights = TRUE,
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 6,
      iter_warmup = 400,
      iter_sampling = 400,
      refresh = 0,
      adapt_delta = 0.9,
      max_treedepth = 12,
      seed = 4102
    )
  )

  expect_true(all(abs(fit$stan_data$marker_weights) < 1e-12))

  assoc_tbl <- summary(fit)$tables$assoc

  alpha_hat <- assoc_tbl$Estimate[assoc_tbl$term == "cv_total"]
  expect_equal(alpha_hat, sim$truth$alpha_cv_total, tolerance = 0.2)

  w_hat <- assoc_tbl$Estimate[grepl("^weight:", assoc_tbl$term)]
  expect_equal(length(w_hat), length(sim$truth$marker_weights))
  expect_equal(w_hat, sim$truth$marker_weights, tolerance = 0.35)
})
