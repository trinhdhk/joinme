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
  pri <- joinme_priors(
    intercept = prior_normal(scale = 2.5),
    slope = prior_student_t(df = 6, scale = 1),
    functional = list(slope = prior_laplace(scale = 0.75)),
    lkj = 2
  )

  expect_s3_class(pri, "joinme_priors")
  expect_equal(pri$global$intercept$scale, 2.5)
  expect_equal(pri$longitudinal$slope$scale, 1.0)
  expect_equal(pri$functional$slope$scale, 0.75)
  expect_equal(pri$lkj$eta, 2)
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

test_that("global roles are inherited independently and formulaDist alpha is unambiguous", {
  pri <- jm_prior(
    intercept = prior_normal(mu = 1, scale = 2),
    slope = prior_student_t(df = 5, mu = 0, scale = 0.8),
    longitudinal = list(slope = prior_laplace(scale = 0.4)),
    alpha = list(intercept = prior_normal(mu = 0, scale = 1.5))
  )

  expect_identical(pri$longitudinal$intercept, pri$global$intercept)
  expect_identical(pri$longitudinal$slope$family, "laplace")
  expect_identical(pri$vcov$sd$slope, pri$global$slope)
  expect_identical(pri$distributional$alpha$intercept$family, "normal")
  expect_identical(pri$distributional$alpha$slope, pri$global$slope)
  expect_error(jm_prior(beta = prior_normal()), "Unknown distributional prior component")
})

test_that("distributional prior selectors normalise quoted family aliases and markers", {
  pri <- jm_prior(
    `sigma[family='student']` = list(
      intercept = prior_laplace(scale = 0.7),
      slope = prior_normal(scale = 0.2)
    ),
    `nu[marker='m2']` = list(intercept = prior_student_t(df = 4))
  )

  expect_identical(
    pri$distributional$sigma$by_family$student_t$intercept$family,
    "laplace"
  )
  expect_equal(
    pri$distributional$sigma$by_family$student_t$slope$scale,
    0.2
  )
  expect_identical(
    pri$distributional$nu$by_marker$m2$intercept$family,
    "student_t"
  )
  expect_equal(pri$distributional$nu$by_marker$m2$intercept$df, 4)
  expect_output(print(pri), "sigma\\[family=student_t\\]\\.slope")
  expect_output(print(pri), "nu\\[marker=m2\\]\\.intercept")
})

test_that("coefficient prior constructors retain separate families and vectors", {
  pri <- jm_prior(
    longitudinal = list(
      intercept = prior_normal(mu = c(0, 1, 2), scale = c(1, 1, 0.5)),
      slope = prior_normal(scale = 0.5)
    ),
    assoc = prior_student_t(df = 4, mu = 0, scale = 1.5),
    functional = list(slope = prior_laplace(mu = 0, scale = 0.75)),
    marker = list(family = prior_horseshoe(df = 3)),
    marker_weights = list(intercept = prior_normal(mu = 0.4, scale = 0.7), family = "normal"),
    vcov = list(
      sd = list(intercept = prior_laplace(mu = c(0, 0.5), scale = c(1, 0.75))),
      corr = list(slope = prior_horseshoe(global_scale = 0.5))
    ),
    class = list(
      slope = prior_laplace(mu = c(0, 1), scale = c(1, 0.5))
    ),
    lkj = prior_lkj(2)
  )

  expect_identical(pri$longitudinal$intercept$family, "normal")
  expect_equal(pri$longitudinal$intercept$mu, c(0, 1, 2))
  expect_identical(pri$assoc$slope$family, "student_t")
  expect_equal(pri$assoc$slope$df, 4)
  expect_identical(pri$functional$slope$family, "laplace")
  expect_identical(pri$marker$family$family, "horseshoe")
  expect_identical(pri$marker_weights$family$family, "normal")
  expect_equal(pri$marker_weights$family$mu, 0)
  expect_equal(pri$marker_weights$family$scale, 1)
  expect_identical(pri$marker_weights$intercept$family, "normal")
  expect_equal(pri$marker_weights$intercept$mu, 0.4)
  expect_equal(pri$marker_weights$intercept$scale, 0.7)
  expect_output(print(pri), "marker_weights.intercept")
  expect_identical(pri$vcov$sd$intercept$family, "laplace")
  expect_equal(pri$vcov$sd$intercept$mu, c(0, 0.5))
  expect_identical(pri$vcov$corr$slope$family, "horseshoe")
  expect_identical(pri$class$slope$family, "laplace")
  expect_equal(pri$class$slope$mu, c(0, 1))
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
    jm_prior(marker_weights = list(family = prior_laplace(mu = 0.5))),
    "family name rather than a prior"
  )
  expect_error(
    jm_prior(marker_weights = list(family = "gaussian")),
    "should be one of"
  )
  marker_horseshoe <- jm_prior(marker = prior_horseshoe(global_scale = 0.2))
  expect_equal(marker_horseshoe$marker$family$global_scale, 0.2)
  expect_error(
    jm_prior(marker_weights = list(family = prior_student_t(df = 8))),
    "family name rather than a prior"
  )
  expect_error(
    jm_prior(marker_weights = prior_student_t(df = 7)),
    "named list"
  )
  marker_weight_intercept <- jm_prior(
    marker_weights = list(intercept = prior_normal(mu = 0.25, scale = 0.6))
  )
  expect_identical(marker_weight_intercept$marker_weights$intercept$family, "normal")
  expect_equal(marker_weight_intercept$marker_weights$intercept$mu, 0.25)
  expect_equal(marker_weight_intercept$marker_weights$intercept$scale, 0.6)
  inherited_marker_weight_intercept <- jm_prior(
    intercept = prior_laplace(mu = -0.1, scale = 0.8),
    marker_weights = list(family = "normal")
  )
  expect_identical(inherited_marker_weight_intercept$marker_weights$intercept$family, "laplace")
  expect_equal(inherited_marker_weight_intercept$marker_weights$intercept$mu, -0.1)
  expect_equal(inherited_marker_weight_intercept$marker_weights$intercept$scale, 0.8)
  expect_error(
    jm_prior(marker_weights = list(mean = prior_normal())),
    "offset.*intercept.*family.*shared"
  )
  named_automatic_marker_weight_student <- jm_prior(
    marker_weights = list(family = "student_t")
  )
  expect_true(named_automatic_marker_weight_student$marker_weights$family$estimate_df)
  expect_output(print(named_automatic_marker_weight_student), "fitted df")
  expect_true(jm_prior()$marker_weights$family$estimate_df)
  expect_error(
    jm_prior(marker_weights = list(family = prior_student_t(mu = 1))),
    "family name rather than a prior"
  )
  expect_error(
    jm_prior(marker_weights = list(family = prior_student_t(scale = 2))),
    "family name rather than a prior"
  )
  expect_error(
    jm_prior(marker_weights = list(family = prior_student_t(df = 2))),
    "family name rather than a prior"
  )
})

test_that("family-only marker blocks do not send ordinary locations or scales to Stan", {
  simulated <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 2),
    times_obs = seq(0, 5, length.out = 2),
    censor_longitudinal_after_event = FALSE,
    seed = 1204,
    use_mirai = FALSE
  )
  stan_data <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = simulated$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = simulated$dataEvent,
    assoc = "cv_total",
    prior_specification = jm_prior(
      marker = list(family = prior_laplace()),
      marker_weights = list(family = "normal")
    )
  )

  expect_false(any(c(
    "prior_marker_mu", "prior_marker_scale",
    "prior_marker_weight_mu", "prior_marker_weight_scale"
  ) %in% names(stan_data)))
  expect_equal(stan_data$prior_marker_family, 3L)
  expect_equal(stan_data$prior_marker_weight_family, 2L)
  expect_equal(stan_data$prior_marker_weight_df, 1)
  expect_equal(stan_data$estimate_marker_weight_df, 0L)
})

test_that("Student-t marker weights fit their degrees of freedom", {
  simulated <- simulate_joinme(
    n_id = 3,
    families = rep("student_t", 3),
    times_obs = seq(0, 5, length.out = 2),
    censor_longitudinal_after_event = FALSE,
    seed = 1205,
    use_mirai = FALSE
  )
  stan_data <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = simulated$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = simulated$dataEvent,
    assoc = "cv_total",
    prior_specification = jm_prior(
      marker_weights = list(
        intercept = prior_laplace(mu = 0.6, scale = 0.3),
        family = "student_t"
      )
    )
  )

  expect_equal(stan_data$prior_marker_weight_family, 1L)
  expect_equal(stan_data$estimate_marker_weight_df, 1L)
  expect_gt(stan_data$n_marker_weight_means, 0L)
  marker_weight_mean_index <- stan_data$prior_start_alpha + 6L + stan_data$M_corr + stan_data$M_vcov
  expect_equal(stan_data$prior_regression_family[marker_weight_mean_index], 3L)
  expect_equal(stan_data$prior_regression_mu[marker_weight_mean_index], 0.6)
  expect_equal(stan_data$prior_regression_scale[marker_weight_mean_index], 0.3)

  automatic_stan_data <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = simulated$dataLong,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    dataEvent = simulated$dataEvent,
    assoc = "cv_total",
    prior_specification = jm_prior(
      marker_weights = list(family = "student_t")
    )
  )
  expect_equal(automatic_stan_data$prior_marker_weight_family, 1L)
  expect_equal(automatic_stan_data$estimate_marker_weight_df, 1L)
})

test_that("marker-weight offset and constant family have one prior interface", {
  named <- jm_prior(marker_weights = list(
    offset = c(m2 = -1, m1 = 2),
    family = "constant"
  ))
  expect_identical(named$marker_weights$family$family, "constant")
  expect_equal(named$marker_weights$offset, c(m2 = -1, m1 = 2))

  synonym <- jm_prior(marker_weights = list(offset = c(2, -1), family = "none"))
  expect_identical(synonym$marker_weights$family$family, "constant")

  expect_true(jm_prior()$marker_weights$shared)
  expect_false(jm_prior(marker_weights = list(shared = FALSE))$marker_weights$shared)
  expect_error(
    jm_prior(marker_weights = list(shared = 0)),
    "marker_weights\\$shared"
  )

  expect_error(
    jm_prior(marker_weights = list(offset = c(m1 = 1, 2))),
    "wholly named or wholly unnamed"
  )
  expect_error(
    jm_prior(marker_weights = list(offset = c(m1 = 1, m1 = 2))),
    "duplicated"
  )
  expect_error(
    jm_prior(marker_weights = list(offset = list(unrecognised = c(1, 2)))),
    "unknown weighted association"
  )
})

test_that("prior assembly is scalar-or-exact rather than partially recycled", {
  pri <- jm_prior(
    longitudinal = list(intercept = prior_normal(mu = c(0, 1), scale = c(1, 2))),
    assoc = prior_student_t(df = 5),
    functional = list(slope = prior_laplace())
  )
  stan_prior <- .pack_regression_priors(
    list(
      beta = list(roles = c("intercept", "intercept"), priors = pri$longitudinal),
      alpha = list(roles = rep("slope", 8L), priors = list(slope = pri$assoc$slope)),
      iota = list(roles = rep("slope", 3L), priors = pri$functional)
    )
  )

  expect_equal(as.integer(stan_prior$prior_regression_family[1:2]), c(2L, 2L))
  expect_equal(stan_prior$prior_regression_mu[1:2], c(0, 1))
  expect_length(stan_prior$prior_regression_scale, 13L)
  expect_equal(as.integer(stan_prior$prior_regression_family[11:13]), rep(3L, 3L))
  expect_error(
    .pack_regression_priors(
      list(beta = list(
        roles = rep("intercept", 3L),
        priors = pri$longitudinal
      ))
    ),
    "exactly one value for each"
  )
})

test_that("horseshoe global scales are shared only within a component-role", {
  pri <- jm_prior(
    longitudinal = list(slope = prior_horseshoe(global_scale = 0.2)),
    vcov = list(sd = list(slope = prior_horseshoe(global_scale = 0.5)))
  )
  packed <- .pack_regression_priors(list(
    longitudinal = list(roles = rep("slope", 2L), priors = pri$longitudinal),
    vcov_sd = list(roles = rep("slope", 3L), priors = pri$vcov$sd)
  ))

  expect_equal(packed$n_regression_horseshoe_local, 5L)
  expect_equal(packed$n_regression_horseshoe_group, 2L)
  expect_equal(as.integer(packed$prior_regression_horseshoe_group_index), c(1L, 1L, 2L, 2L, 2L))
  expect_equal(packed$prior_regression_horseshoe_global_scale, c(0.2, 0.5))
})
