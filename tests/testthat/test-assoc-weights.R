test_that("association summaries include marker weights when estimated", {
  set.seed(123)
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 4,
    D = 2,
    n_t = 3,
    seed = 123,
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
  expect_true(any(grepl("^weight:", assoc_tbl$term)))

  re <- ranef(fit)
  expect_true(!is.null(re$assoc_weight))
  expect_true(any(grepl("^weight:", re$assoc_weight$term)))
})
