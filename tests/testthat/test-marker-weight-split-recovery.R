test_that("fit supports end-to-end split marker-weight structures", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  term_weights <- list(
    cv_total = c(1.0, 0.4),
    cv_marker = c(0.2, 1.2)
  )
  assoc_truth <- c(cv_total = 0.55, cv_marker = 0.35)

  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    n_id = 50,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 5,
    times_obs = seq(0, 7, length.out = 10),
    marker_weights = term_weights,
    shared_marker_weights = FALSE,
    assoc = c("cv_total", "cv_marker"),
    assoc_coefs = assoc_truth,
    transforms = list(
      cv_total = list(type = "identity"),
      cv_marker = list(type = "identity")
    ),
    seed = 6021
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total", "cv_marker"),
    families = rep("gaussian", 2),
    marker_weights = term_weights,
    shared_marker_weights = FALSE,
    fixed_marker_weights = TRUE,
    transforms = list(
      cv_total = list(type = "identity"),
      cv_marker = list(type = "identity")
    ),
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 2,
      iter_warmup = 120,
      iter_sampling = 120,
      refresh = 0,
      adapt_delta = 0.9,
      max_treedepth = 12,
      seed = 6021
    )
  )

  assoc_tbl <- summary(fit)$tables$assoc
  cv_row <- assoc_tbl[assoc_tbl$term == "cv_total (+)", , drop = FALSE]
  cs_row <- assoc_tbl[assoc_tbl$term == "cv_marker (+)", , drop = FALSE]

  expect_equal(nrow(cv_row), 1L)
  expect_equal(nrow(cs_row), 1L)
  expect_true(is.finite(cv_row$Estimate) && cv_row$Estimate > 0)
  expect_true(is.finite(cs_row$Estimate) && cs_row$Estimate > 0)

  assoc_draws <- assoc(fit, summary = FALSE)
  expect_true(is.matrix(assoc_draws$cv_total))
  expect_true(is.matrix(assoc_draws$cv_marker))

  assoc_map <- extract(fit, what = "assoc", keep_chains = FALSE)$term_map
  expect_true(any(grepl("^weight\\[cv_total\\]:", assoc_map$term)))
  expect_true(any(grepl("^weight\\[cv_marker\\]:", assoc_map$term)))
})
