test_that("renamed public class aliases preserve old S3 compatibility", {
  fit <- JoiNMe::JoiNMeFit$new(
    fit = NULL,
    stan_data = list(),
    formulaLong = y ~ 1 + time,
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = NULL,
    config = list(),
    call = quote(JoiNMe::joinme(formulaLong = y ~ 1 + time)),
    tmax = 1,
    dataLong = data.frame(id = c(1, 1), time = c(0, 1), marker = "m1", y = c(1, 2)),
    dataEvent = data.frame(id = 1, time = 1.2, event = 0L)
  )

  pred <- JoiNMe::PredJoiNMeFit$new(
    predictions = list(),
    quantiles = list(),
    draws = list(),
    data = list(),
    metadata = list(),
    call = quote(predict(fit)),
    tmax = 1,
    n_samples = 1
  )

  expect_true(inherits(fit, "JoiNMeFit"))
  expect_true(inherits(fit, "JoiNMeFit"))
  expect_true(inherits(pred, "PredJoiNMeFit"))
  expect_true(inherits(pred, "JoiNMeDynPred"))
})
