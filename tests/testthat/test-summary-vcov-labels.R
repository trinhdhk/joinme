test_that("print.summary_JoinMeFit uses bespoke covariance labels", {
  s <- SummaryJoinMeFit$new(
    tables = list(
      fixef = NULL,
      gamma_w = NULL,
      assoc = NULL,
      distributional = NULL,
      distributional_regression = NULL,
      vcov = list(
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
