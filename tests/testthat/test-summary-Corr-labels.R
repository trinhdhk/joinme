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
        marker = data.frame(block = "marker", row = 1L, col = 1L, Estimate = 2)
      ),
      id_marker_cov = list(
        latent = data.frame(
          block = "sigma_latent",
          row = c(1L, 2L),
          col = c(1L, 1L),
          Estimate = c(0.7, 0.1),
          Est.Error = c(0.1, 0.05),
          Q2.5 = c(0.5, -0.1),
          Q97.5 = c(0.9, 0.25),
          Rhat = rep(1.0, 2),
          ess_bulk = rep(100, 2),
          ess_tail = rep(100, 2)
        ),
        regression = data.frame(
          block = c("L", "L", "L", "global"),
          row = c(1L, 1L, 1L, NA_integer_),
          col = c(1L, 1L, 1L, NA_integer_),
          term = c("(Intercept)", "x1", "lambda", "sd_u"),
          Estimate = c(-0.2, 0.05, 0.3, 0.6),
          Est.Error = c(0.1, 0.03, 0.1, 0.08),
          Q2.5 = c(-0.4, -0.02, 0.1, 0.45),
          Q97.5 = c(0.0, 0.12, 0.5, 0.75),
          Rhat = rep(1.0, 4),
          ess_bulk = rep(100, 4),
          ess_tail = rep(100, 4)
        )
      )
    ),
    diagnostics = NULL,
    metadata = list(family = "student_t", tmax = 1)
  )

  txt <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(txt, "\\nid\\n")
  expect_match(txt, "\\nmarker\\n")
  expect_true(grepl("id:marker covariance parameters", txt, fixed = TRUE))
  expect_true(grepl("latent covariance matrix", txt, fixed = TRUE))
  expect_true(grepl("covariance regression coefficients", txt, fixed = TRUE))
  expect_true(grepl("sigma_latent", txt, fixed = TRUE))
  expect_true(grepl("block", txt, fixed = TRUE))
  expect_true(grepl("row", txt, fixed = TRUE))
  expect_true(grepl("col", txt, fixed = TRUE))
  expect_true(grepl("\\(Intercept\\)", txt))
  expect_true(grepl("x1", txt, fixed = TRUE))
  expect_true(grepl("lambda", txt, fixed = TRUE))
  expect_true(grepl("sd_u", txt, fixed = TRUE))
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
