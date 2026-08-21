test_that(".get_stan_file always resolves the threaded Stan sources", {
  expect_warning(
    fit_file <- .get_stan_file("joinme_fit", threaded = FALSE),
    "deprecated"
  )
  pred_file <- .get_stan_file("joinme_dynpred", threaded = TRUE)
  mix_fit_file <- .get_stan_file(
    "joinme_mix_fit",
    threaded = TRUE
  )
  mix_pred_file <- .get_stan_file(
    "joinme_mix_dynpred",
    threaded = TRUE
  )
  fitpred_file <- .get_stan_file("joinme_fitpred", threaded = TRUE)
  mix_fitpred_file <- .get_stan_file(
    "joinme_mix_fitpred",
    threaded = TRUE
  )

  expect_match(basename(fit_file), "joinme_fit_threading\\.stan$")
  expect_match(basename(pred_file), "joinme_dynpred_threading\\.stan$")
  expect_match(basename(mix_fit_file), "joinme_mix_fit_threading\\.stan$")
  expect_match(
    basename(mix_pred_file),
    "joinme_mix_dynpred_threading\\.stan$"
  )
  expect_match(basename(fitpred_file), "joinme_fitpred_threading\\.stan$")
  expect_match(
    basename(mix_fitpred_file),
    "joinme_mix_fitpred_threading\\.stan$"
  )
})

test_that("Stan programme routing follows the fitted model family", {
  expect_identical(
    .stan_fit_program(list(use_mixture = 0L)),
    "joinme_fit"
  )
  expect_identical(
    .stan_fit_program(list(use_mixture = 1L)),
    "joinme_mix_fit"
  )

  ordinary <- list(stan_data = list(use_mixture = 0L))
  mixture <- structure(
    list(stan_data = list(use_mixture = 1L)),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )
  expect_identical(
    .stan_dynpred_program(ordinary),
    "joinme_dynpred"
  )
  expect_identical(
    .stan_dynpred_program(mixture),
    "joinme_mix_dynpred"
  )
  expect_identical(
    .stan_fitpred_program(ordinary),
    "joinme_fitpred"
  )
  expect_identical(
    .stan_fitpred_program(mixture),
    "joinme_mix_fitpred"
  )
})

test_that(".get_stan_engine falls back to an available backend", {
  testthat::local_mocked_bindings(
    .stan_backend_available = function(engine) identical(engine, "rstan"),
    .package = "joinme"
  )

  expect_warning(
    resolved <- .get_stan_engine("cmdstanr"),
    "Falling back"
  )
  expect_identical(resolved, "rstan")
})

test_that("single-line Stan declarations explain their role inline", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  stan_files <- list.files(
    stan_root,
    pattern = "\\.(stan|stanfunctions)$",
    recursive = TRUE,
    full.names = TRUE
  )
  declaration_pattern <- paste0(
    "^\\s*(array\\[[^;]+\\]\\s+)?",
    "(int|real|vector|row_vector|matrix|simplex|ordered|",
    "positive_ordered|unit_vector|cholesky_factor_corr|",
    "corr_matrix|cov_matrix)",
    "(<[^;]+>)?(\\[[^;]+\\])?\\s+",
    "[A-Za-z][A-Za-z0-9_]*(\\s*=\\s*[^;]+)?;\\s*$"
  )
  unexplained <- character(0)
  for (stan_file in stan_files) {
    source_lines <- readLines(stan_file, warn = FALSE)
    declaration_lines <- grep(
      declaration_pattern,
      source_lines,
      perl = TRUE
    )
    missing_comment <- declaration_lines[
      !grepl("//", source_lines[declaration_lines], fixed = TRUE)
    ]
    if (length(missing_comment) > 0L) {
      unexplained <- c(
        unexplained,
        paste0(
          sub(paste0("^", stan_root, "/"), "", stan_file),
          ":",
          missing_comment
        )
      )
    }
  }
  expect_identical(unexplained, character(0))
})
