test_that("standard normal CDF and probit quantile have distinct bytecodes", {
  normal_ops <- joinme:::.bytecode_normal_ops()
  expect_identical(normal_ops, c(PHI = 26L, INV_PHI = 27L))

  cdf_programs <- list(
    Phi = parse_transform_expr(~ Phi(x)),
    pnorm = parse_transform_expr(~ pnorm(x))
  )
  quantile_programs <- list(
    inv_Phi = parse_transform_expr(~ inv_Phi(x)),
    qnorm = parse_transform_expr(~ qnorm(x)),
    probit = parse_transform_expr(~ probit(x))
  )

  for (program in cdf_programs) {
    expect_identical(program$bytecode, c(0L, normal_ops[["PHI"]]))
    expect_equal(
      joinme:::eval_bytecode_vector(c(-2, 0, 2), program$bytecode),
      stats::pnorm(c(-2, 0, 2))
    )
  }

  probabilities <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  for (program in quantile_programs) {
    expect_identical(program$bytecode, c(0L, normal_ops[["INV_PHI"]]))
    expect_equal(
      joinme:::eval_bytecode_vector(probabilities, program$bytecode),
      stats::qnorm(probabilities)
    )
  }
})

test_that("formula-link inversion respects the two normal directions", {
  probit_formula <- invert_transform_expr(~ inv_Phi(x))
  probit_alias_formula <- invert_transform_expr(~ probit(x))
  cdf_formula <- invert_transform_expr(~ Phi(x))

  expect_equal(
    joinme:::eval_bytecode_vector(
      c(-1, 0, 1),
      parse_transform_expr(probit_formula)$bytecode
    ),
    stats::pnorm(c(-1, 0, 1))
  )
  expect_equal(
    joinme:::eval_bytecode_vector(
      c(-1, 0, 1),
      parse_transform_expr(probit_alias_formula)$bytecode
    ),
    stats::pnorm(c(-1, 0, 1))
  )

  probabilities <- c(0.1, 0.5, 0.9)
  expect_equal(
    joinme:::eval_bytecode_vector(
      probabilities,
      parse_transform_expr(cdf_formula)$bytecode
    ),
    stats::qnorm(probabilities)
  )
})

test_that("normal bytecodes retain fitted affine-shift indices", {
  transform <- joinme_tf(
    cv_total = ~ Phi(x, intercept = TRUE, slope = TRUE)
  )
  program <- parse_transform_expr(
    transform$cv_total$expr,
    iota_nodes = transform$cv_total$iota_nodes
  )

  expect_identical(program$bytecode, c(0L, 26L))
  expect_identical(program$op_iota_intercept_idx, c(0L, 1L))
  expect_identical(program$op_iota_slope_idx, c(0L, 1L))
  expect_equal(
    joinme:::eval_bytecode_vector(
      c(-1, 0, 1),
      bytecode = program$bytecode,
      iota_intercepts = 0.4,
      iota_slopes = 1.5,
      op_iota_intercept_idx = program$op_iota_intercept_idx,
      op_iota_slope_idx = program$op_iota_slope_idx
    ),
    stats::pnorm(0.4 + 1.5 * c(-1, 0, 1))
  )
})

test_that("named and formula probit families compile the Phi inverse link", {
  named <- jm_family("bernoulli", link = "probit")
  formula <- jm_family("bernoulli", link = ~ inv_Phi(x))
  direct <- jm_family("bernoulli", inv_link = ~ Phi(x))

  for (specification in list(named, formula, direct)) {
    expect_identical(specification$link, "probit")
    expect_identical(specification$inv_link$bytecode, c(0L, 26L))
    expect_equal(
      joinme:::eval_bytecode_vector(
        c(-1, 0, 1),
        specification$inv_link$bytecode
      ),
      stats::pnorm(c(-1, 0, 1))
    )
  }
})

test_that("response-scale methods share one inverse-link implementation", {
  eta <- matrix(c(-1, 0, 1, 2), nrow = 2)

  expect_equal(
    joinme:::.apply_inverse_link_matrix(eta, 1L),
    eta
  )
  expect_equal(
    joinme:::.apply_inverse_link_matrix(eta, 2L),
    exp(eta)
  )
  expect_equal(
    joinme:::.apply_inverse_link_matrix(eta, 3L),
    matrix(stats::plogis(eta), nrow = nrow(eta))
  )
  expect_equal(
    joinme:::.apply_inverse_link_matrix(eta, 4L),
    matrix(stats::pnorm(eta), nrow = nrow(eta))
  )

  positive_eta <- exp(eta)
  expect_equal(
    joinme:::.apply_inverse_link_matrix(positive_eta, 5L),
    eta
  )

  probabilities <- matrix(c(0.1, 0.25, 0.75, 0.9), nrow = 2)
  expect_equal(
    joinme:::.apply_inverse_link_matrix(
      probabilities,
      link_code = 0L,
      bytecode = c(0L, 27L)
    ),
    matrix(stats::qnorm(probabilities), nrow = nrow(probabilities))
  )
})

test_that("simulation evaluates the inverse-normal bytecode instruction", {
  # qnorm(plogis(eta)) is finite for every finite eta. This custom inverse link
  # therefore exercises instruction 27 in the simulation evaluator without
  # imposing an invalid probability-domain assumption on the linear predictor.
  custom_family <- jm_family(
    "gaussian",
    inv_link = ~ inv_Phi(inv_logit(x))
  )

  simulation <- simulate_joinme(
    n_id = 4,
    families = rep(list(custom_family), 3),
    n_obs_per_marker_per_id = 2,
    times_obs = c(0, 1),
    time_cens = 1.5,
    seed = 2701,
    use_mirai = FALSE
  )

  expect_true(nrow(simulation$dataLong) > 0L)
  expect_true(all(is.finite(simulation$dataLong$y)))
  expect_true(27L %in% as.integer(custom_family$inv_link$bytecode))
})

test_that("fit and dynamic-prediction Stan data admit both normal operations", {
  root <- normalizePath(
    if (file.exists(file.path("inst", "stan", "helper", "data", "fit_data.stan"))) {
      "."
    } else {
      file.path("..", "..")
    },
    mustWork = TRUE
  )
  stan_files <- file.path(
    root,
    "inst", "stan", "helper", "data",
    c("fit_data.stan", "dynpred_data.stan")
  )
  stan_text <- lapply(stan_files, readLines, warn = FALSE)

  expect_true(all(vapply(
    stan_text,
    function(lines) any(grepl("inv_link_ops.*upper=27", lines, fixed = FALSE)),
    logical(1)
  )))

  evaluator <- readLines(
    file.path(
      root, "inst", "stan", "helper", "functions",
      "functional_transform.stanfunctions"
    ),
    warn = FALSE
  )
  expect_true(any(grepl("op == 26", evaluator, fixed = TRUE)))
  expect_true(any(grepl("Phi(", evaluator, fixed = TRUE)))
  expect_true(any(grepl("op == 27", evaluator, fixed = TRUE)))
  expect_true(any(grepl("inv_Phi(", evaluator, fixed = TRUE)))
})
