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
    times_obs = seq(0, 8, length.out = 12),
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.6))),
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

  expect_true(all(abs(fit$stan_data$marker_weight_offsets) < 1e-12))

  assoc_tbl <- summary(fit)$tables$assoc

  alpha_hat <- assoc_tbl$Estimate[assoc_tbl$term == "cv_total (+)"]
  expect_true(is.finite(alpha_hat) && alpha_hat > 0)

  w_rows <- marker_weights(fit)
  truth_w <- sim$truth$marker_weights
  marker_names <- sim$marker_info$names
  if (is.null(names(truth_w))) names(truth_w) <- marker_names
  w_terms <- as.character(w_rows$marker)
  expect_equal(length(w_terms), length(truth_w))
  for (i in seq_along(w_terms)) {
    w_true <- truth_w[[w_terms[i]]]
    expect_true(!is.na(w_true))
    expect_true(w_true >= w_rows$Q2.5[i] && w_true <= w_rows$Q97.5[i])
  }
})
