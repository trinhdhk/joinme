test_that("Stan entry points assemble scientific submodels", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  include_root <- file.path(stan_root, "include")
  submodel_root <- file.path(include_root, "submodels")
  components <- c(
    "longitudinal", "survival", "assoc", "marker_weight",
    "functional", "latent_class"
  )

  expect_true(file.exists(file.path(stan_root, "README.md")))
  expect_true(file.exists(file.path(include_root, "README.md")))
  expect_true(file.exists(file.path(submodel_root, "README.md")))
  expect_true(all(file.exists(file.path(submodel_root, components, "README.md"))))

  stan_blocks <- c(
    "data", "transformed_data", "parameters", "transformed_parameters",
    "model", "generated_quantities"
  )
  block_paths <- unlist(lapply(
    components,
    function(component) file.path(submodel_root, component, stan_blocks)
  ))
  expect_true(all(dir.exists(block_paths)))
  expect_true(all(file.exists(file.path(block_paths, "README.md"))))

  entry_points <- c(
    "joinme_fit_threading.stan", "joinme_mix_fit_threading.stan",
    "joinme_dynpred_threading.stan", "joinme_mix_dynpred_threading.stan"
  )
  fitpred_entry_points <- c(
    "joinme_fitpred_threading.stan",
    "joinme_mix_fitpred_threading.stan"
  )
  entry_text <- lapply(file.path(stan_root, entry_points), readLines, warn = FALSE)
  fitpred_text <- lapply(
    file.path(stan_root, fitpred_entry_points),
    readLines,
    warn = FALSE
  )
  expect_true(all(file.exists(file.path(stan_root, fitpred_entry_points))))
  expect_true(all(vapply(fitpred_text, function(lines) {
    any(grepl("include/etc/data/fitted_random_effect_draw_data.stan", lines, fixed = TRUE))
  }, logical(1))))
  expect_false(dir.exists(file.path(submodel_root, "full_model")))
  common_components <- setdiff(components, "latent_class")
  expect_true(all(vapply(entry_text, function(lines) {
    all(vapply(
      common_components,
      function(component) any(grepl(
        paste0("include/submodels/", component, "/"),
        lines,
        fixed = TRUE
      )),
      logical(1)
    ))
  }, logical(1))))
  expect_false(any(grepl("include/submodels/latent_class/", entry_text[[1L]], fixed = TRUE)))
  expect_true(any(grepl("include/submodels/latent_class/", entry_text[[2L]], fixed = TRUE)))
  expect_false(any(grepl("include/submodels/latent_class/", entry_text[[3L]], fixed = TRUE)))
  expect_true(any(grepl("include/submodels/latent_class/", entry_text[[4L]], fixed = TRUE)))
  expect_false(dir.exists(file.path(stan_root, "helper")))
  expect_false(dir.exists(file.path(stan_root, "submodels")))
})

test_that("every submodel Stan fragment starts with a Doxygen file account", {
  submodel_root <- testthat::test_path(
    "..", "..", "inst", "stan", "include", "submodels"
  )
  stan_files <- list.files(
    submodel_root,
    pattern = "[.]stan$",
    recursive = TRUE,
    full.names = TRUE
  )

  expect_gt(length(stan_files), 0L)
  has_file_account <- vapply(stan_files, function(path) {
    initial_lines <- readLines(path, n = 24L, warn = FALSE)
    any(grepl("@file", initial_lines, fixed = TRUE)) &&
      any(grepl("@brief", initial_lines, fixed = TRUE)) &&
      any(grepl("@details", initial_lines, fixed = TRUE))
  }, logical(1))

  expect_true(all(has_file_account), info = paste(stan_files[!has_file_account], collapse = "\n"))
})

test_that("cross-submodel Stan code lives under include/etc", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  helper_paths <- c(
    file.path("include", "etc", "data", "regression_prior_data.stan"),
    file.path("include", "etc", "data", "dynamic_prediction_draw_data.stan"),
    file.path("include", "etc", "parameters", "regression_prior_parameters.stan"),
    file.path("include", "etc", "parameters", "dynamic_prediction_subject_effects.stan"),
    file.path("include", "etc", "transformed_parameters", "regression_coefficients.stan"),
    file.path("include", "etc", "model", "regression_priors.stan")
  )
  expect_true(all(file.exists(file.path(stan_root, helper_paths))))
})

test_that("Stan data discovery ignores semicolons in block explanations", {
  stan_file <- testthat::test_path("..", "..", "inst", "stan", "joinme_fit_threading.stan")
  data_names <- .stan_data_names(stan_file)

  expect_true("estimate_marker_weights" %in% data_names)
  expect_false("zero" %in% data_names)
  expect_false("case" %in% data_names)
})
