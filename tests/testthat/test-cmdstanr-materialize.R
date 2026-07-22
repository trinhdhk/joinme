test_that(".import_cmdstanr_fit eagerly loads cmdstanr-backed contents", {
  calls <- character(0)

  fake_fit <- structure(list(
    draws = function(...) {
      calls <<- c(calls, "draws")
      invisible(matrix(1, 1, 1))
    },
    sampler_diagnostics = function(...) {
      calls <<- c(calls, "sampler_diagnostics")
      stop("diagnostics unavailable")
    },
    init = function(...) {
      calls <<- c(calls, "init")
      invisible(list())
    },
    profiles = function(...) {
      calls <<- c(calls, "profiles")
      invisible(list())
    }
  ), class = "CmdStanMCMC")

  out <- joinme:::.import_cmdstanr_fit(fake_fit)

  expect_identical(out, fake_fit)
  expect_equal(calls, c("draws", "sampler_diagnostics", "init", "profiles"))
})

test_that(".import_cmdstanr_fit leaves non-cmdstan fits unchanged", {
  fit_like <- list(a = 1)
  out <- joinme:::.import_cmdstanr_fit(fit_like)
  expect_identical(out, fit_like)
})
