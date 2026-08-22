test_that("named links compile their mathematical inverse links", {
  log_spec <- jm_family("poisson", link = "log")
  exp_spec <- jm_family("gaussian", link = "exp")

  expect_equal(log_spec$link, "log")
  expect_equal(log_spec$inv_link$bytecode, c(0L, 7L))
  expect_equal(.link_name_from_inv_link_expr(~ exp(x)), "log")
  expect_equal(
    eval_bytecode_vector(
      c(-1, 0, 1),
      log_spec$inv_link$bytecode,
      log_spec$inv_link$const_data
    ),
    exp(c(-1, 0, 1))
  )

  expect_equal(exp_spec$link, "exp")
  expect_equal(exp_spec$inv_link$bytecode, c(0L, 6L))
  expect_equal(.link_name_from_inv_link_expr(~ log(x)), "exp")
  expect_equal(
    eval_bytecode_vector(
      c(1, exp(1)),
      exp_spec$inv_link$bytecode,
      exp_spec$inv_link$const_data
    ),
    c(0, 1)
  )
})

test_that("jm_family prints inverse links as mathematical expressions", {
  named_output <- capture.output(print(jm_family("poisson", link = "log")))
  expect_true(any(grepl("Family: poisson", named_output, fixed = TRUE)))
  expect_true(any(grepl("Link: log", named_output, fixed = TRUE)))
  expect_true(any(grepl("Inverse link: exp(x)", named_output, fixed = TRUE)))

  custom_output <- capture.output(print(jm_family(
    "gaussian",
    inv_link = ~ inv_logit(2 * x - 0.25)
  )))
  expect_true(any(grepl("Link: custom", custom_output, fixed = TRUE)))
  expect_true(any(grepl(
    "Inverse link: inv_logit(2 * x - 0.25)",
    custom_output,
    fixed = TRUE
  )))

  formula_output <- capture.output(print(jm_family(
    "gaussian",
    link = ~ log(2 * x + 1)
  )))
  expect_true(any(grepl(
    "Inverse link: (exp(x) - 1)/2",
    formula_output,
    fixed = TRUE
  )))
})

test_that("inverse-link expression reconstruction names every unary instruction", {
  unary_labels <- c(
    `6` = "log(x)",
    `7` = "exp(x)",
    `8` = "sqrt(x)",
    `9` = "inv_logit(x)",
    `10` = "logit(x)",
    `11` = "1/x",
    `13` = "sin(x)",
    `14` = "cos(x)",
    `15` = "tan(x)",
    `16` = "abs(x)",
    `17` = "x^2",
    `18` = "sinh(x)",
    `19` = "cosh(x)",
    `20` = "tanh(x)",
    `21` = "asinh(x)",
    `22` = "acosh(x)",
    `23` = "atanh(x)",
    `24` = "log1p_exp(x)",
    `25` = "cbrt(x)",
    `26` = "Phi(x)",
    `27` = "inv_Phi(x)"
  )

  reconstructed_labels <- vapply(names(unary_labels), function(operation) {
    base::deparse1(.bytecode_as_expression(c(0L, as.integer(operation))))
  }, character(1))
  expect_identical(unname(reconstructed_labels), unname(unary_labels))
})

test_that("formula links are inverted before bytecode compilation", {
  simple <- jm_family("poisson", link = ~ log(x))
  affine <- jm_family("gaussian", link = ~ log(2 * x + 1))
  cubic <- jm_family("gaussian", link = ~ x^3)

  expect_equal(simple$link, "log")
  expect_equal(simple$inv_link$bytecode, c(0L, 7L))
  expect_equal(
    eval_bytecode_vector(
      log(c(3, 5)),
      affine$inv_link$bytecode,
      affine$inv_link$const_data
    ),
    c(1, 2)
  )
  expect_equal(
    eval_bytecode_vector(
      c(-8, -1, 0, 1, 8, 27),
      cubic$inv_link$bytecode,
      cubic$inv_link$const_data
    ),
    c(-2, -1, 0, 1, 2, 3)
  )

  expect_error(
    jm_family("gaussian", link = ~ x + x),
    "exactly once"
  )
  expect_error(
    jm_family("gaussian", link = ~ abs(x)),
    "not invertible"
  )
  expect_error(
    jm_family("gaussian", link = ~ x^2),
    "not one-to-one"
  )
})

test_that("formula links pass through standata for fit and dynamic prediction", {
  set.seed(191)
  sim <- simulate_joinme(
    n_id = 8,
    families = c("gaussian", "poisson"),
    times_obs = seq(0, 3, length.out = 5),
    assoc = "cv_total",
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2))),
    seed = 191
  )
  standata <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = "cv_total",
    families = list(
      jm_family("gaussian", link = ~ x),
      jm_family("poisson", link = ~ log(x))
    )
  )

  expect_equal(standata$link_names, c("identity", "log"))
  expect_equal(as.integer(standata$link_long), c(1L, 2L))
  expect_equal(as.integer(standata$inv_link_ops[2, 1:2]), c(0L, 7L))

  # The fitted object retains these arrays and predict.JoiNMeFit copies them
  # into dynpred standata, so the same inverse-link program is used by both
  # Stan programs.
  expect_equal(as.integer(standata$inv_link_n_ops), c(1L, 2L))
  expect_equal(as.numeric(standata$inv_link_const[2, ]), 0)
})
