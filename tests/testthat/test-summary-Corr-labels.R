test_that("print.summary_JoinMeFit uses bespoke covariance labels", {
  s <- SummaryJoinMeFit$new(
    tables = list(
      fixef = NULL,
      gamma_w = NULL,
      assoc = NULL,
      distributional = NULL,
      distributional_regression = NULL,
      corr = list(
        id = data.frame(block = "id", row = 1L, col = 1L, Estimate = 1),
        marker = data.frame(block = "marker", row = 1L, col = 1L, Estimate = 2),
        marker_by_id_latent = data.frame(block = "id:marker", row = 1L, col = 1L, Estimate = 3)
      )
    ),
    diagnostics = NULL,
    metadata = list(family = "student_t", tmax = 1)
  )

  txt <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(txt, "\\nid\\n")
  expect_match(txt, "\\nmarker\\n")
  expect_match(txt, "\\nid:marker\\n")
  expect_false(grepl("Sigma_u|Sigma_v|Sigma_w", txt))
})

test_that("print.summary_JoinMeFit prints survival process report only when available", {
  s_with <- SummaryJoinMeFit$new(
    tables = list(
      fixef = NULL,
      gamma_w = NULL,
      survival_process = data.frame(
        term = "x1",
        Estimate = 0.2,
        Hazard.Ratio = 1.22,
        Est.Error = 0.1,
        Q2.5 = -0.1,
        Q97.5 = 0.5,
        HR.Q2.5 = 0.90,
        HR.Q97.5 = 1.65,
        Rhat = 1.00,
        ess_bulk = 100,
        ess_tail = 100,
        stringsAsFactors = FALSE
      ),
      assoc = NULL,
      distributional = NULL,
      distributional_regression = NULL,
      corr = NULL
    ),
    diagnostics = NULL,
    metadata = list(family = "student_t", tmax = 1)
  )

  txt_with <- paste(capture.output(print(s_with)), collapse = "\n")
  expect_true(grepl("Survival process (non-association covariates)", txt_with, fixed = TRUE))

  s_without <- SummaryJoinMeFit$new(
    tables = list(
      fixef = NULL,
      gamma_w = NULL,
      survival_process = NULL,
      assoc = NULL,
      distributional = NULL,
      distributional_regression = NULL,
      corr = NULL
    ),
    diagnostics = NULL,
    metadata = list(family = "student_t", tmax = 1)
  )

  txt_without <- paste(capture.output(print(s_without)), collapse = "\n")
  expect_false(grepl("Survival process (non-association covariates)", txt_without, fixed = TRUE))
})
