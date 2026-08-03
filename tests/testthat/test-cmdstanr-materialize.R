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

test_that("chain energy diagnostics use the retained sampler order", {
  sampler_array <- array(
    NA_real_,
    dim = c(4L, 2L, 6L),
    dimnames = list(
      iteration = as.character(1:4),
      chain = as.character(1:2),
      variable = c(
        "treedepth__",
        "divergent__",
        "energy__",
        "accept_stat__",
        "stepsize__",
        "n_leapfrog__"
      )
    )
  )
  sampler_array[, , "treedepth__"] <- 3
  sampler_array[, , "divergent__"] <- 0
  sampler_array[2L, 1L, "divergent__"] <- 1
  sampler_array[, 1L, "energy__"] <- c(1, 2, 1, 3)
  sampler_array[, 2L, "energy__"] <- c(2, 4, 3, 5)
  sampler_array[, , "accept_stat__"] <- 0.9
  sampler_array[, , "stepsize__"] <- 0.1
  sampler_array[, , "n_leapfrog__"] <- 7
  fake_fit <- structure(
    list(
      sampler_diagnostics = function(format = "draws_array", ...) {
        posterior::as_draws_array(sampler_array)
      }
    ),
    class = "CmdStanMCMC"
  )

  energy <- joinme:::.sampler_energy_by_chain(fake_fit)

  expect_equal(nrow(energy), 2L)
  expect_equal(energy$divergences, c(1L, 0L))
  expect_equal(
    energy$E_BFMI[1L],
    mean(diff(c(1, 2, 1, 3))^2) / stats::var(c(1, 2, 1, 3))
  )
  expect_equal(energy$median_leapfrog, c(7, 7))
  expect_equal(energy$step_size, c(0.1, 0.1))
})
