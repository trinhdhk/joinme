test_that("standardized shrinkage draws match all three Stan families", {
  set.seed(7301)
  got_t <- .sim_draw_standard_shrinkage(8, 0L)
  set.seed(7301)
  expect_equal(got_t, stats::rt(8, df = 6))

  set.seed(7302)
  got_laplace <- .sim_draw_standard_shrinkage(8, 1L)
  set.seed(7302)
  expected_laplace <- sample(c(-1, 1), 8, replace = TRUE) * stats::rexp(8, rate = 1)
  expect_equal(got_laplace, expected_laplace)

  set.seed(7303)
  got_normal <- .sim_draw_standard_shrinkage(8, 2L)
  set.seed(7303)
  expect_equal(got_normal, stats::rnorm(8))
})

test_that("simulate_joinme stores base, mean, departure, and effective marker-weight truth", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 3),
    times_obs = c(0, 0.5, 1),
    time_cens = 1,
    shrinkage = 1L,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0)),
      marker_weights = list(
        offset = c(m1 = 0.5, m2 = -0.25, m3 = 1),
        intercept = prior_normal(mu = 0.2, scale = 0.4),
        family = "laplace"
      )
    ),
    assoc = "cv_total",
    seed = 7304,
    use_mirai = FALSE
  )

  truth <- sim$truth
  expect_false("true_params" %in% names(sim))
  expect_equal(truth$shrinkage, 1L)
  expect_equal(truth$shrinkage_distribution, "double_exponential(0, 1)")
  expect_equal(truth$marker_weight_prior_family, "laplace")
  expect_identical(truth$marker_weight_prior_family, "laplace")
  expect_equal(
    unname(truth$marker_weights),
    unname(truth$marker_weights_offset + truth$marker_weight_mean[[1L]] + truth$marker_weights_latent)
  )
  expect_equal(
    unname(truth$stan_fit$z_marker_weights),
    unname(truth$marker_weights_latent)
  )
  expect_null(truth$marker_weight_scale)
  expect_equal(
    unname(truth$stan_fit$marker_weight_standardised),
    unname(truth$marker_weights_latent)
  )
  expect_true(any(abs(truth$marker_weights_latent) > 0))
  expect_equal(
    truth$recovery$arguments$priors$marker_weights$offset,
    c(m1 = 0.5, m2 = -0.25, m3 = 1)
  )
  expect_equal(truth$recovery$arguments$priors$marker_weights$intercept$mu, 0.2)
  expect_equal(truth$recovery$arguments$priors$marker_weights$intercept$scale, 0.4)
  expect_equal(truth$recovery$arguments$priors$marker_weights$family$mu, 0)
  expect_equal(truth$recovery$arguments$priors$marker_weights$family$scale, 1)
})

test_that("marker-weight population means are fixed when supplied and drawn once when omitted", {
  simulation_arguments <- list(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.2),
    time_cens = 0.3,
    truth = jm_truth(assoc_coef = list(slope = c(cv_marker = 0)),
      marker_weights = list(family = "normal")
    ),
    assoc = "cv_marker",
    use_mirai = FALSE
  )

  drawn_a <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7321)))
  drawn_b <- do.call(simulate_joinme, c(simulation_arguments, list(seed = 7322)))
  fixed_arguments <- simulation_arguments
  fixed_arguments$truth <- jm_truth(
    assoc_coef = list(slope = c(cv_marker = 0)),
    marker_weights = list(intercept = 1.25, family = "normal")
  )
  fixed <- do.call(simulate_joinme, c(fixed_arguments, list(seed = 7323)))

  expect_length(drawn_a$truth$marker_weight_mean, 1L)
  expect_true(is.finite(drawn_a$truth$marker_weight_mean[["shared"]]))
  expect_false(identical(
    drawn_a$truth$marker_weight_mean,
    drawn_b$truth$marker_weight_mean
  ))
  expect_identical(fixed$truth$marker_weight_mean, c(shared = 1.25))
  expect_equal(fixed$truth$marker_weight_mean_by_term$cv_marker, 1.25)
})

test_that("simulate_joinme uses zero offsets for estimated weights and constant offsets exactly", {
  estimated <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    shrinkage = 2L,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0)),
      marker_weights = list(family = "normal")
    ),
    assoc = "cv_total",
    seed = 7305,
    use_mirai = FALSE
  )
  expect_equal(unname(estimated$truth$marker_weights_offset), c(0, 0))
  expect_equal(estimated$truth$shrinkage_distribution, "normal(0, 1)")
  expect_equal(estimated$truth$marker_weight_prior_family, "normal")

  supplied <- c(m1 = 2, m2 = -1)
  fixed <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0)),
      marker_weights = list(offset = supplied, family = "constant")
    ),
    shrinkage = 0L,
    assoc = "cv_total",
    seed = 7306,
    use_mirai = FALSE
  )
  expect_equal(fixed$truth$marker_weights, supplied)
  expect_identical(fixed$truth$marker_weight_prior_family, "constant")
  expect_equal(unname(fixed$truth$marker_weights_latent), c(0, 0))
  expect_equal(fixed$truth$stan_fit$z_marker_weights, numeric(0))
})

test_that("term-specific simulated weights use Stan's flattened set layout", {
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0, cv_marker = 0)),
      marker_weights = list(
        offset = list(cv_total = c(0.2, 0.4), cv_marker = c(-0.3, 0.1)),
        intercept = c(cv_total = 0.7, cv_marker = -0.4),
        shared = FALSE
      )
    ),
    shrinkage = 0L,
    assoc = c("cv_total", "cv_marker"),
    seed = 7307,
    use_mirai = FALSE
  )

  truth <- sim$truth
  expect_equal(names(truth$marker_weights), c("cv_total", "cv_marker"))
  expect_equal(
    truth$stan_fit$z_marker_weights,
    as.numeric(t(truth$marker_weight_standardised_sets))
  )
  for (term in c("cv_total", "cv_marker")) {
    expect_equal(
      unname(truth$marker_weights[[term]]),
      unname(truth$marker_weights_offset[[term]] + truth$marker_weight_mean_by_term[[term]] + truth$marker_weights_latent[[term]])
    )
  }
  expect_equal(truth$marker_weight_mean, c(cv_total = 0.7, cv_marker = -0.4))
  expect_length(truth$marker_weight_df, 1L)
  expect_gt(truth$marker_weight_df, 2)
  expect_identical(truth$marker_weight_df_was_fitted, TRUE)
  expect_null(truth$marker_weight_scale)
})

test_that("mixture simulation retains the same direct marker-weight departures", {
  simulation <- simulate_joinme_mix(
    n_id = 4,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.2),
    time_cens = 0.3,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0)),
      marker_weights = list(
        offset = c(m1 = 0.4, m2 = -0.2),
        family = "normal"
      )
    ),
    assoc = "cv_total",
    n_classes = 2,
    class_type = "subject",
    class_dimensions = 1,
    seed = 7308,
    use_mirai = FALSE
  )

  truth <- simulation$truth
  expect_equal(truth$z_marker_weight_sets, truth$marker_weight_standardised_sets)
  expect_equal(
    truth$stan_fit$z_marker_weights,
    as.numeric(t(truth$z_marker_weight_sets))
  )
  expect_true(any(abs(truth$z_marker_weight_sets) > 0))
  expect_false("marker_weight_scale" %in% names(truth$stan_fit))
  expect_false("true_params" %in% names(simulation))
  expect_length(truth$marker_weight_mean, 1L)
  expect_true(is.finite(truth$marker_weight_mean[["shared"]]))
  expect_equal(
    truth$marker_weight_mean_by_term$cv_total,
    truth$marker_weight_mean[["shared"]]
  )
})

test_that("mixture truth retains term-specific marker-weight means", {
  simulation <- simulate_joinme_mix(
    n_id = 4,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.2),
    time_cens = 0.3,
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0, cs_marker = 0)),
      marker_weights = list(
        intercept = c(cv_total = 0.65, cs_marker = -0.35),
        shared = FALSE,
        family = "normal"
      )
    ),
    assoc = c("cv_total", "cs_marker"),
    n_classes = 2,
    class_type = "subject",
    class_dimensions = 1,
    seed = 7309,
    use_mirai = FALSE
  )

  expect_equal(
    simulation$truth$marker_weight_mean,
    c(cv_total = 0.65, cs_marker = -0.35)
  )
  expect_equal(
    unlist(
      simulation$truth$marker_weight_mean_by_term[c("cv_total", "cs_marker")],
      use.names = TRUE
    ),
    c(cv_total = 0.65, cs_marker = -0.35)
  )
  expect_false("true_params" %in% names(simulation))
})

test_that("simulate_joinme validates marker-weight declarations", {
  expect_error(
    jm_prior(marker_weights = list(offset = c(m1 = 1, 2))),
    "wholly named or wholly unnamed"
  )
  expect_error(
    simulate_joinme(n_id = 2, shrinkage = 3L),
    "shrinkage"
  )
  expect_error(
    simulate_joinme(n_id = 2, marker_weight_scale = 1, use_mirai = FALSE),
    "unused argument"
  )
  expect_error(
    simulate_joinme_mix(n_id = 2, marker_weight_scale = 1, use_mirai = FALSE),
    "unused argument"
  )
  fixed <- simulate_joinme(
    n_id = 2,
    truth = jm_truth(marker_weights = list(
      offset = c(m1 = 2, m2 = -1, m3 = 0.5),
      family = "none"
    )),
    seed = 7310,
    use_mirai = FALSE
  )
  expect_equal(fixed$truth$marker_weights, c(m1 = 2, m2 = -1, m3 = 0.5))
})

test_that("Stan uses direct unit-scale marker-weight departures", {
  parameter_code <- paste(
    readLines(testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "marker_weight", "parameters", "fit.stan"), warn = FALSE),
    collapse = "\n"
  )
  transformed_code <- paste(
    readLines(testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "marker_weight", "transformed_parameters", "fit.stan"), warn = FALSE),
    collapse = "\n"
  )
  model_code <- paste(
    readLines(testthat::test_path("..", "..", "inst", "stan", "include", "submodels", "marker_weight", "model", "fit.stan"), warn = FALSE),
    collapse = "\n"
  )

  expect_false(grepl("marker_weight_scale", parameter_code, fixed = TRUE))
  expect_false(grepl("marker_weight_scale", model_code, fixed = TRUE))
  expect_match(
    parameter_code,
    "n_marker_weight_means > 0 ? 1 : 0] marker_weight_df_excess",
    fixed = TRUE
  )
  expect_match(
    model_code,
    "z_marker_weights ~ student_t(marker_weight_df[1], 0, 1);",
    fixed = TRUE
  )
  expect_false(grepl("marker_weight_df[s]", model_code, fixed = TRUE))
  expect_match(
    transformed_code,
    "z_marker_weight_sets[s] = to_row_vector(marker_weight_prior_effect[start_pos:end_pos]);",
    fixed = TRUE
  )
  expect_false(grepl("marker_weight_scale", transformed_code, fixed = TRUE))
})

test_that("regularised-horseshoe marker-weight simulation retains its realised hierarchy", {
  set.seed(7311)
  draw <- .sim_draw_marker_weight_prior(
    5L,
    prior_horseshoe(global_scale = 0.25, slab_scale = 1)
  )

  expect_length(draw$departure, 5L)
  expect_length(draw$raw, 5L)
  expect_length(draw$local_scale, 5L)
  expect_length(draw$global_scale, 1L)
  expect_length(draw$slab_multiplier, 1L)
  expect_true(all(is.finite(draw$departure)))
})

test_that("Student-t marker-weight simulation draws fitted degrees of freedom", {
  set.seed(7312)
  moving_prior <- jm_prior(marker_weights = list(family = "student_t"))$marker_weights$family
  automatic_draw <- .sim_draw_marker_weight_prior(
    12L,
    moving_prior,
    n_sets = 3L
  )

  expect_length(automatic_draw$departure, 12L)
  expect_length(automatic_draw$df, 1L)
  expect_true(all(automatic_draw$df > 2))
  expect_identical(automatic_draw$df_was_fitted, TRUE)
})
