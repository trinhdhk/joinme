test_that("stan cache key tracks include dependency content", {
  tmp <- tempfile("joinme-stan-key-")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  root_file <- file.path(tmp, "model.stan")
  helper_dir <- file.path(tmp, "helper")
  dir.create(helper_dir, recursive = TRUE, showWarnings = FALSE)
  helper_file <- file.path(helper_dir, "f.stanfunctions")

  writeLines(c(
    "functions {",
    "  #include helper/f.stanfunctions",
    "}",
    "data { int<lower=0> N; }",
    "parameters { real y; }",
    "model { y ~ normal(0, 1); }"
  ), root_file)

  writeLines(c(
    "real f1(real x) {",
    "  return x;",
    "}"
  ), helper_file)

  key1 <- joinme:::.stan_cache_key(root_file, cpp_options = list(stan_threads = TRUE))

  writeLines(c(
    "real f1(real x) {",
    "  return x + 1e-9;",
    "}"
  ), helper_file)

  key2 <- joinme:::.stan_cache_key(root_file, cpp_options = list(stan_threads = TRUE))

  expect_true(nzchar(key1))
  expect_true(nzchar(key2))
  expect_false(identical(key1, key2))
})

test_that("stan cache key reflects cpp options", {
  tmp <- tempfile("joinme-stan-key-cpp-")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  root_file <- file.path(tmp, "model_simple.stan")

  writeLines(c(
    "data { int<lower=0> N; }",
    "parameters { real y; }",
    "model { y ~ normal(0, 1); }"
  ), root_file)

  key_no_threads <- joinme:::.stan_cache_key(root_file, cpp_options = list(stan_threads = FALSE))
  key_threads <- joinme:::.stan_cache_key(root_file, cpp_options = list(stan_threads = TRUE))

  expect_false(identical(key_no_threads, key_threads))
})

test_that("stan cache key is stable across different absolute root paths", {
  tmp1 <- tempfile("joinme-stan-key-path-a-")
  tmp2 <- tempfile("joinme-stan-key-path-b-")
  dir.create(tmp1, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp2, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp1, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(tmp2, recursive = TRUE, force = TRUE), add = TRUE)

  make_tree <- function(root) {
    root_file <- file.path(root, "model.stan")
    helper_dir <- file.path(root, "helper")
    dir.create(helper_dir, recursive = TRUE, showWarnings = FALSE)
    helper_file <- file.path(helper_dir, "f.stanfunctions")

    writeLines(c(
      "functions {",
      "  #include helper/f.stanfunctions",
      "}",
      "data { int<lower=0> N; }",
      "parameters { real y; }",
      "model { y ~ normal(0, 1); }"
    ), root_file)

    writeLines(c(
      "real f1(real x) {",
      "  return x;",
      "}"
    ), helper_file)

    root_file
  }

  root1 <- make_tree(tmp1)
  root2 <- make_tree(tmp2)

  key1 <- joinme:::.stan_cache_key(root1, cpp_options = list(stan_threads = TRUE))
  key2 <- joinme:::.stan_cache_key(root2, cpp_options = list(stan_threads = TRUE))

  expect_equal(key1, key2)
})
