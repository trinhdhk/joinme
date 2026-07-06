test_that(".get_stan_file always resolves the threaded Stan sources", {
  expect_warning(
    fit_file <- JoiNMe:::.get_stan_file("JoiNMe_fit", threaded = FALSE),
    "deprecated"
  )
  pred_file <- JoiNMe:::.get_stan_file("JoiNMe_dynpred", threaded = TRUE)

  expect_match(basename(fit_file), "JoiNMe_fit_threading\\.stan$")
  expect_match(basename(pred_file), "JoiNMe_dynpred_threading\\.stan$")
})

test_that(".resolve_stan_engine falls back to an available backend", {
  testthat::local_mocked_bindings(
    .stan_backend_available = function(engine) identical(engine, "rstan"),
    .package = "joinme"
  )

  expect_warning(
    resolved <- JoiNMe:::.resolve_stan_engine("cmdstanr"),
    "Falling back"
  )
  expect_identical(resolved, "rstan")
})
