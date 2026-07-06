testthat::test_that("print.JoiNMeDynPred shows key metadata", {
  pred <- JoiNMe::JoiNMeDynPred$new(
    predictions = list(longitudinal = NULL, survival = NULL, cumhaz = NULL),
    quantiles = list(),
    draws = list(),
    data = list(),
    metadata = list(pred_type = "per_marker_id", scale = "epred", n_subjects = 2),
    call = quote(predict(fit_obj)),
    tmax = 12,
    n_samples = 200
  )

  out <- paste(capture.output(print(pred)), collapse = "\n")

  testthat::expect_match(out, "JoiNMe dynamic prediction", fixed = TRUE)
  testthat::expect_match(out, "Prediction type: per_marker_id", fixed = TRUE)
  testthat::expect_match(out, "Scale: epred", fixed = TRUE)
  testthat::expect_match(out, "Subjects: 2", fixed = TRUE)
  testthat::expect_match(out, "Posterior draws: 200", fixed = TRUE)
  testthat::expect_match(out, "tmax: 12", fixed = TRUE)
})
