test_that("joinme_tf validates and normalises shorthand declarations", {
  tf <- joinme_tf(
    cv_total = "identity",
    corr = ~ -x,
    cv_marker = list(type = "pwlin", x = c(-1, 0, 1), y = c(0.2, 1, 0.2))
  )

  expect_s3_class(tf, "joinme_tf")
  expect_identical(tf$cv_total$type, "identity")
  expect_identical(tf$corr$type, "functional")
  expect_true(inherits(tf$corr$expr, "formula"))
  expect_identical(tf$cv_marker$type, "pwlin")

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$tf_mode_cv_tot, 0)
  expect_identical(sd_tf$tf_mode_corr, 1)
  expect_identical(sd_tf$tf_mode_cv_marker, 3)
})

test_that("joinme_tf supports expit-based spline transforms", {
  expit_grid <- stats::plogis(seq(-4, 4, length.out = 21))
  tf <- joinme_tf(
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      n_knots = 5,
      degree = 2,
      lambda = 0.5
    )
  )

  expect_s3_class(tf, "joinme_tf")
  expect_identical(tf$corr$type, "ispline_expit_penalised")

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$tf_mode_corr, 4L)
  expect_true(all(sd_tf$knots_corr > 0 & sd_tf$knots_corr < 1))
  expect_equal(length(sd_tf$coeff_corr), sd_tf$n_coeff_corr)
})

test_that("joinme_tf rejects raw-scale spline inputs for expit-based transforms", {
  expect_error(
    joinme_tf(
      corr = list(
        type = "ispline_expit_penalised",
        x = seq(-4, 4, length.out = 21),
        n_knots = 5,
        degree = 2,
        lambda = 0.5
      )
    ),
    "must lie on the expit scale"
  )
})

test_that("joinme_tf rejects unknown channels", {
  expect_error(
    joinme_tf(not_a_channel = "identity"),
    "Unknown transform channel"
  )
})

test_that("joinme_tf records fit-only affine shift flags for supported functional transforms", {
  tf <- joinme_tf(cv_total = ~ SoftMax(x, intercept = TRUE, slope = TRUE))

  expect_identical(tf$cv_total$type, "functional")
  expect_identical(as.character(tf$cv_total$expr[[2]][[1]]), "softmax")
  expect_true(isTRUE(tf$cv_total$estimate_iota_intercept))
  expect_true(isTRUE(tf$cv_total$estimate_iota_slope))

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$estimate_iota_intercept_cv, 1L)
  expect_identical(sd_tf$estimate_iota_slope_cv, 1L)
})

test_that("joinme_tf tracks multiple affine-shifted nonlinear nodes", {
  tf <- joinme_tf(
    cv_total = ~ softplus(expit(x, intercept = TRUE, slope = TRUE), intercept = TRUE, slope = TRUE)
  )

  expect_identical(tf$cv_total$n_iota_intercept, 2L)
  expect_identical(tf$cv_total$n_iota_slope, 2L)
  expect_length(tf$cv_total$iota_nodes, 2L)
  expect_identical(vapply(tf$cv_total$iota_nodes, `[[`, character(1), "path"), c("root/1", "root"))

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$estimate_iota_intercept_cv, 2L)
  expect_identical(sd_tf$estimate_iota_slope_cv, 2L)
  expect_equal(sd_tf$functional_iota_intercept_idx_cv, c(0L, 1L, 2L))
  expect_equal(sd_tf$functional_iota_slope_idx_cv, c(0L, 1L, 2L))
})

test_that("joinme_tf rejects affine shift flags outside supported functional transforms", {
  expect_error(
    joinme_tf(cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 0.5, 1), intercept = TRUE)),
    "Fit-only affine-shift options"
  )
})

test_that("joinme_priors validates exposed prior components", {
  pri <- joinme_priors(beta = list(scale = 2.5), alpha = 1.0, iota = 0.75, lkj = 2)

  expect_s3_class(pri, "joinme_priors")
  expect_equal(pri$beta$scale, 2.5)
  expect_equal(pri$alpha$scale, 1.0)
  expect_equal(pri$iota$scale, 0.75)
  expect_equal(pri$lkj$eta, 2)

  stan_pri <- .build_priors(beta_prior = pri$beta, alpha_prior = pri$alpha, iota_prior = pri$iota, lkj_prior = pri$lkj)
  expect_equal(stan_pri$alpha$scale, 1.0)
  expect_equal(stan_pri$iota$scale, 0.75)
  expect_equal(stan_pri$lkj_eta, 2)
})

test_that("joinme_priors rejects invalid LKJ values", {
  expect_error(
    joinme_priors(lkj = 0),
    "positive finite"
  )
})

test_that("singular prior aliases share the validated entry point", {
  expect_identical(jm_prior, joinme_priors)
  expect_identical(joinme_prior, joinme_priors)
})

test_that("coefficient prior constructors retain separate families and vectors", {
  pri <- jm_prior(
    beta = prior_normal(mu = c(0, 1, 2), scale = c(1, 1, 0.5)),
    alpha = prior_student_t(df = 4, mu = 0, scale = 1.5),
    iota = prior_laplace(mu = 0, scale = 0.75),
    marker = prior_horseshoe(df = 3),
    marker_weight = prior_normal(),
    vcov_sd = prior_laplace(mu = c(0, 0.5), scale = c(1, 0.75)),
    vcov_corr = prior_horseshoe(global_scale = 0.5),
    class_regression = prior_laplace(mu = c(0, 1), scale = c(1, 0.5)),
    lkj = prior_lkj(2)
  )

  expect_identical(pri$beta$family, "normal")
  expect_equal(pri$beta$mu, c(0, 1, 2))
  expect_identical(pri$alpha$family, "student_t")
  expect_equal(pri$alpha$df, 4)
  expect_identical(pri$iota$family, "laplace")
  expect_identical(pri$marker$family, "horseshoe")
  expect_identical(pri$marker_weight$family, "normal")
  expect_identical(pri$vcov_sd$family, "laplace")
  expect_equal(pri$vcov_sd$mu, c(0, 0.5))
  expect_identical(pri$vcov_corr$family, "horseshoe")
  expect_identical(pri$class_regression$family, "laplace")
  expect_equal(pri$class_regression$mu, c(0, 1))
  expect_equal(pri$lkj$eta, 2)
})

test_that("marker_random_effect is no longer a prior component", {
  expect_false("marker_random_effect" %in% names(formals(joinme_priors)))
  expect_error(
    .joinme_priors_(list(marker_random_effect = prior_normal())),
    "Unknown prior component"
  )
})

test_that("coefficient prior constructors accept scalar numeric lists", {
  declared_prior <- prior_normal(
    mu = list(0, 1, 2, 3),
    scale = list(1, 1, 2, 0.5)
  )

  expect_equal(declared_prior$mu, c(0, 1, 2, 3))
  expect_equal(declared_prior$scale, c(1, 1, 2, 0.5))
  expect_error(
    prior_normal(mu = list(c(0, 1), 2)),
    "list of single numeric values"
  )
})

test_that("coefficient prior degrees of freedom are fixed scalars", {
  expect_error(
    prior_student_t(df = c(3, 6)),
    "must be one positive finite number"
  )
  expect_error(
    prior_horseshoe(global_df = c(1, 2)),
    "must be one positive finite number"
  )
})

test_that("standardised marker priors reject user-supplied location and scale", {
  expect_error(
    jm_prior(marker = prior_normal(scale = 2)),
    "standardised latent block"
  )
  expect_error(
    jm_prior(marker_weight = prior_laplace(mu = 0.5)),
    "standardised latent block"
  )
  expect_error(
    jm_prior(marker = prior_horseshoe(global_scale = 0.2)),
    "standardised latent block"
  )
})

test_that("prior materialisation is scalar-or-exact rather than partially recycled", {
  encoded <- .build_priors(
    beta_prior = prior_normal(mu = c(0, 1), scale = c(1, 2)),
    alpha_prior = prior_student_t(df = 5),
    iota_prior = prior_laplace(),
    marker_prior = prior_normal(),
    marker_weight_prior = prior_normal(),
    lkj_prior = prior_lkj(1.5)
  )
  stan_prior <- .materialise_joinme_prior_data(
    encoded,
    c(beta = 2L, alpha = 8L, iota = 3L, marker = 4L, marker_weight = 2L)
  )

  expect_identical(stan_prior$prior_beta_family, 2L)
  expect_equal(stan_prior$prior_beta_mu, c(0, 1))
  expect_length(stan_prior$prior_alpha_scale, 8L)
  expect_identical(stan_prior$prior_iota_family, 3L)
  expect_error(
    .materialise_joinme_prior_data(
      encoded,
      c(beta = 3L, alpha = 8L, iota = 3L, marker = 4L, marker_weight = 2L)
    ),
    "exactly one value per parameter"
  )
})
