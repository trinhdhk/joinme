test_that("marker-weight summaries are compact and individual weights have an accessor", {
  set.seed(123)
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("student_t", 2),
    times_obs = seq(0, 4, length.out = 8),
    seed = 123,
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.6))),
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
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      refresh = 0,
      adapt_delta = 0.95,
      max_treedepth = 12,
      seed = 123
    )
  )

  s <- summary(fit)
  assoc_tbl <- s$tables$assoc
  expect_false(any(grepl("^weight:", assoc_tbl$term)))
  expect_equal(s$tables$marker_weights$term, c("mean weight", "SD weight"))

  individual <- marker_weights(fit)
  expect_equal(individual$marker, fit$stan_data$marker_levels)
  expect_true(all(c("offset", "Estimate", "Est.Error", "Q2.5", "Q97.5") %in% names(individual)))

  re <- ranef(fit)
  expect_null(re$formulaLong$marker_weight)
  expect_true(!is.null(re$assoc))
  expect_true(all(c("assoc_term", "marker") %in% names(re$assoc)))
})
