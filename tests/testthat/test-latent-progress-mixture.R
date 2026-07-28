test_that("ordinary models receive a parameter-free mixture representation", {
  mixture_data <- .build_mixture_standata(
    stan_data = list(),
    mixture = NULL
  )

  expect_identical(mixture_data$use_mixture, 0L)
  expect_identical(mixture_data$n_clusters, 1L)
  expect_identical(mixture_data$K_mix, 0L)
  expect_length(mixture_data$mix_idx_subject, 0L)
  expect_identical(mixture_data$P_class_subject, 0L)
  expect_identical(mixture_data$P_class_marker, 0L)
  expect_null(mixture_data$mixture)
})

test_that("compatible random-effect levels share G rather than forming a product", {
  stan_data <- list(
    R_id = 3L,
    R_mk = 2L,
    Q_idm = 2L,
    indep_idmarker_cov = 0L,
    estimate_marker_weights = 0L,
    use_marker_weight_assoc = 0L,
    n_marker_weight_sets = 0L,
    shrinkage = 0L,
    shared_marker_weights = 1L
  )
  mixture_data <- .build_mixture_standata(
    stan_data,
    mixture = list(
      n_clusters = 3L,
      cluster_type = c("subject", "vcov"),
      cluster_dimensions = NULL,
      class_probability_concentration = 1,
      class_regression_scale = 1,
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$n_clusters, 3L)
  expect_identical(mixture_data$K_mix, 4L)
  expect_identical(mixture_data$mix_idx_subject, 1:2)
  expect_identical(mixture_data$mix_idx_covariance, 1:2)
  expect_identical(mixture_data$mix_start_subject, 1L)
  expect_identical(mixture_data$mix_start_covariance, 3L)
  expect_identical(
    mixture_data$mixture$allocation_domains$subject,
    c("subject", "vcov")
  )
  expect_equal(mixture_data$mix_probability_prior, rep(1, 3))
})

test_that("corr selects only off-diagonal covariance-regression coordinates", {
  mixture_data <- .build_mixture_standata(
    stan_data = list(
      R_id = 1L,
      R_mk = 1L,
      Q_idm = 3L,
      indep_idmarker_cov = 0L,
      shrinkage = 2L
    ),
    mixture = list(
      n_clusters = 2L,
      cluster_type = "corr",
      cluster_dimensions = c(1L, 2L),
      class_probability_concentration = 1,
      class_regression_scale = 1,
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$mix_idx_covariance, c(2L, 4L))
  expect_identical(mixture_data$mixture$cluster_type, "corr")
  expect_identical(mixture_data$mixture$dimensions$corr, c(2L, 4L))
  expect_false("mix_marker_weight" %in% names(mixture_data))
})

test_that("subject and marker domains retain common labels and separate units", {
  stan_data <- list(
    R_id = 2L,
    R_mk = 2L,
    Q_idm = 1L,
    indep_idmarker_cov = 1L,
    estimate_marker_weights = 0L,
    use_marker_weight_assoc = 0L,
    n_marker_weight_sets = 0L,
    shrinkage = 2L,
    shared_marker_weights = 1L
  )
  mixture_data <- .build_mixture_standata(
    stan_data,
    mixture = list(
      n_clusters = 2L,
      cluster_type = c("subject", "marker"),
      cluster_dimensions = NULL,
      class_probability_concentration = c(2, 3),
      class_regression_scale = 1,
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$n_clusters, 2L)
  expect_identical(mixture_data$K_mix, 4L)
  expect_equal(mixture_data$mix_probability_prior, c(2, 3))
  expect_identical(
    mixture_data$mixture$allocation_domains,
    list(subject = "subject", marker = "marker")
  )
})

test_that("cluster_type accepts only the four public types", {
  expect_identical(
    .canonical_cluster_types(c("subject", "marker", "corr")),
    c("subject", "marker", "corr")
  )
  expect_error(.canonical_cluster_types("individual"), "subject")
  expect_error(.canonical_cluster_types("id_marker"), "subject")
  expect_error(.canonical_cluster_types("marker_weight"), "subject")
  expect_error(
    .canonical_cluster_types(c("corr", "vcov")),
    "cannot contain both"
  )
  expect_error(
    joinme_mix(
      y ~ time,
      data.frame(id = 1L, marker = "m", time = 0, y = 0),
      cluster = "subject",
      fit = FALSE
    ),
    "matches multiple formal arguments"
  )
})

test_that("mixture dimensions and unavailable levels are explained", {
  stan_data <- list(
    R_id = 1L,
    R_mk = 0L,
    Q_idm = 0L,
    indep_idmarker_cov = 1L,
    estimate_marker_weights = 0L,
    use_marker_weight_assoc = 0L,
    n_marker_weight_sets = 0L,
    shrinkage = 0L
  )
  specification <- list(
    n_clusters = 2L,
    cluster_type = "subject",
    cluster_dimensions = 2L,
    class_probability_concentration = 1,
    class_regression_scale = 1,
    include_survival = FALSE
  )

  expect_error(
    .build_mixture_standata(stan_data, specification),
    "Choose distinct integer indices between 1 and 1"
  )
  specification$cluster_type <- "marker"
  specification$cluster_dimensions <- NULL
  expect_error(
    .build_mixture_standata(stan_data, specification),
    "unavailable"
  )
  specification$cluster_type <- "subject"
  specification$n_clusters <- 2.5
  expect_error(
    .build_mixture_standata(stan_data, specification),
    "single integer"
  )
})

test_that("longitudinal-only scaffolds preserve identifiers and follow-up", {
  longitudinal_data <- data.frame(
    id = c("b", "a", "b", "a"),
    time = c(0, 0, 2, 1),
    marker = c("m1", "m1", "m1", "m1"),
    y = 1:4
  )
  scaffold <- .mixture_event_scaffold(
    longitudinal_data,
    id_variable = "id",
    time_variable = "time"
  )

  expect_equal(scaffold$id, c("a", "b"))
  expect_equal(scaffold$.joinme_mix_follow_up, c(2, 2))
  expect_true(all(scaffold$.joinme_mix_event == 0L))
})

test_that("formulaCluster creates unit-level designs and stable probabilities", {
  longitudinal_data <- data.frame(
    id = rep(1:3, each = 2L),
    marker = rep(c("m1", "m2"), times = 3L),
    treatment = rep(c(0, 1, 0), each = 2L),
    marker_type = factor(
      rep(c("reference", "active"), times = 3L),
      levels = c("reference", "active")
    )
  )
  event_data <- data.frame(
    id = 1:3,
    treatment = c(0, 1, 0)
  )
  design <- .mixture_class_design(
    formula_class = list(
      subject = ~ treatment,
      marker = ~ marker_type
    ),
    selected_types = c("subject", "marker"),
    data_long = longitudinal_data,
    data_event = event_data,
    id_variable = "id",
    marker_variable = "marker"
  )

  expect_identical(design$subject$columns, "class_1:treatment")
  expect_identical(design$marker$columns, "class_1:marker_typeactive")
  expect_equal(design$subject$matrix[, 1L], c(0, 1, 0))
  expect_equal(design$marker$matrix[, 1L], c(0, 1))

  baseline <- matrix(c(0.4, 0.6, 0.4, 0.6), nrow = 2L, byrow = TRUE)
  coefficients <- matrix(c(log(2), log(3)), nrow = 2L, ncol = 1L)
  probability <- .mixture_class_probability(
    baseline_probability = baseline,
    class_coefficient = coefficients,
    class_design = matrix(c(0, 1), ncol = 1L),
    class_term_start = c(1L, 1L),
    class_term_count = c(1L, 0L)
  )
  expect_equal(probability[, 1L, ], baseline)
  expect_true(all(probability[, 2L, 1L] > baseline[, 1L]))
  expect_equal(apply(probability, c(1L, 2L), sum), matrix(1, 2L, 2L))
})

test_that("class-specific formulae use compact coefficients and no default ordering", {
  longitudinal_data <- data.frame(
    id = rep(1:3, each = 2L),
    marker = "m1",
    age = rep(c(20, 30, 40), each = 2L),
    treatment = rep(c(0, 1, 0), each = 2L)
  )
  event_data <- data.frame(
    id = 1:3,
    age = c(20, 30, 40),
    treatment = c(0, 1, 0)
  )
  design <- .mixture_class_design(
    formula_class = list(~age, ~treatment, ~age + treatment),
    selected_types = "subject",
    data_long = longitudinal_data,
    data_event = event_data,
    id_variable = "id",
    marker_variable = "marker",
    n_clusters = 3L
  )

  expect_true(design$class_specific)
  expect_identical(
    design$subject$class_term_count,
    c(1L, 1L, 2L)
  )
  expect_identical(
    design$subject$class_term_start,
    c(1L, 2L, 3L)
  )
  expect_identical(
    design$subject$coefficient_class,
    c(1L, 2L, 3L, 3L)
  )
  expect_identical(
    .resolve_cluster_ordering(NULL, design$class_specific),
    "none"
  )
  expect_identical(
    .resolve_cluster_ordering("probability", design$class_specific),
    "probability"
  )
})

test_that("cluster_ordering targets only an intercept, baseline probabilities, or neither", {
  stan_data <- list(
    n_id = 3L,
    D = 1L,
    R_id = 2L,
    R_mk = 0L,
    Q_idm = 0L,
    zid_cols = c("(Intercept)", "time"),
    zmk_cols = character(0),
    zidm_cols = character(0),
    indep_idmarker_cov = 1L,
    estimate_marker_weights = 0L,
    use_marker_weight_assoc = 0L,
    n_marker_weight_sets = 0L,
    shrinkage = 0L
  )
  specification <- list(
    n_clusters = 3L,
    cluster_type = "subject",
    cluster_dimensions = c(1L, 2L),
    class_probability_concentration = 1,
    class_regression_scale = 1,
    cluster_ordering = "intercept",
    include_survival = FALSE
  )

  intercept_ordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(intercept_ordered$mix_ordering, 1L)
  expect_identical(
    intercept_ordered$mix_ordered_location_coordinate,
    1L
  )

  specification$cluster_ordering <- "probability"
  probability_ordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(probability_ordered$mix_ordering, 2L)
  expect_identical(
    probability_ordered$mix_ordered_location_coordinate,
    0L
  )

  specification$cluster_ordering <- "none"
  unordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(unordered$mix_ordering, 0L)
  expect_identical(unordered$mix_ordered_location_coordinate, 0L)
})

test_that("latent-class priors are declared through jm_prior", {
  priors <- jm_prior(
    class_probability = c(2, 3, 4),
    class_regression = 1.25
  )
  expect_s3_class(priors, "joinme_priors")
  expect_equal(priors$class_probability, c(2, 3, 4))
  expect_identical(priors$class_regression, 1.25)
  expect_false("class_probability_prior" %in% names(formals(joinme_mix)))
  expect_false("class_location_scale" %in% names(formals(joinme_mix)))
  expect_false("class_scale_rate" %in% names(formals(joinme_mix)))
  expect_false("class_separation" %in% names(formals(joinme_mix)))
})

test_that("joinme_mix prepares the inherited holder through the common entry path", {
  captured_arguments <- NULL
  fake_holder <- new.env(parent = emptyenv())
  fake_holder$stan_data <- list(mixture = list(
    n_clusters = 2L,
    cluster_type = "subject",
    include_survival = TRUE
  ))
  fake_holder$sample <- function() "sampled mixture"

  testthat::local_mocked_bindings(
    joinme = function(...) {
      captured_arguments <<- list(...)
      fake_holder
    },
    .package = "joinme"
  )

  prepared <- joinme_mix(
    y ~ time,
    data.frame(id = 1L, time = 0, marker = "m", y = 1),
    survival::Surv(time, event) ~ 1,
    data.frame(id = 1L, time = 1, event = 0L),
    n_clusters = 2,
    cluster_type = "subject",
    fit = FALSE
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(prepared$result_class, "mixture")
  expect_identical(captured_arguments$fit, FALSE)
  expect_identical(captured_arguments$mixture$n_clusters, 2L)
  expect_identical(captured_arguments$mixture$cluster_type, "subject")
  expect_identical(captured_arguments$mixture$cluster_ordering, "intercept")
})

test_that("joinme_mix defaults class-specific formulae to no ordering", {
  captured_arguments <- NULL
  fake_holder <- new.env(parent = emptyenv())
  fake_holder$stan_data <- list(mixture = list(
    n_clusters = 2L,
    levels = "subject",
    include_survival = TRUE
  ))
  fake_holder$sample <- function() "sampled mixture"
  testthat::local_mocked_bindings(
    joinme = function(...) {
      captured_arguments <<- list(...)
      fake_holder
    },
    .package = "joinme"
  )

  prepared <- joinme_mix(
    y ~ time,
    data.frame(
      id = c(1L, 2L),
      time = 0,
      marker = "m",
      x = c(0, 1),
      y = c(1, 2)
    ),
    survival::Surv(time, event) ~ 1,
    data.frame(
      id = c(1L, 2L),
      time = 1,
      event = 0L,
      x = c(0, 1)
    ),
    n_clusters = 2,
    cluster_type = "subject",
    formulaCluster = list(~x, ~x),
    fit = FALSE
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(captured_arguments$mixture$cluster_ordering, "none")
  expect_identical(
    captured_arguments$mixture$class_design$subject$class_term_count,
    c(1L, 1L)
  )
})

test_that("Stan fit data and source include every mixture field", {
  stan_file <- file.path(
    system.file(package = "joinme"),
    "stan",
    "joinme_mix_fit_threading.stan"
  )
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stan_file <- testthat::test_path(
      "..",
      "..",
      "inst",
      "stan",
      "joinme_mix_fit_threading.stan"
    )
  }
  required_data <- .stan_data_names(stan_file)
  expect_true(all(c(
    "use_mixture",
    "n_clusters",
    "K_mix",
    "mix_idx_subject",
    "mix_idx_marker",
    "mix_idx_covariance",
    "mix_ordering",
    "mix_ordered_location_coordinate",
    "mix_probability_prior",
    "P_class_subject",
    "X_class_subject",
    "class_term_start_subject",
    "class_term_count_subject",
    "P_class_marker",
    "X_class_marker",
    "class_term_start_marker",
    "class_term_count_marker",
    "class_regression_scale"
  ) %in% required_data))

  source <- paste(.read_stan_with_includes(stan_file), collapse = "\n")
  expect_match(source, "fit_mixture_priors", fixed = TRUE)
  expect_match(source, "ordered[n_clusters] mix_location_ordered", fixed = TRUE)
  expect_match(source, "mix_location_unordered", fixed = TRUE)
  expect_match(source, "positive_ordered[n_clusters - 1]", fixed = TRUE)
  expect_match(source, "posterior_class_probability_subject", fixed = TRUE)
  expect_match(source, "posterior_class_probability_marker", fixed = TRUE)
})

test_that("ordinary Stan fitting has no latent-class contract", {
  stan_file <- testthat::test_path(
    "..",
    "..",
    "inst",
    "stan",
    "joinme_fit_threading.stan"
  )
  required_data <- .stan_data_names(stan_file)
  source <- paste(.read_stan_with_includes(stan_file), collapse = "\n")

  expect_false(any(grepl("^mix_|^n_clusters$|^K_mix$|^use_mixture$", required_data)))
  expect_false(grepl("latent_progress_component_lpdf", source, fixed = TRUE))
  expect_false(grepl("posterior_class_probability_", source, fixed = TRUE))
})

test_that("dynamic prediction restores the fitted mixture distribution", {
  stan_file <- file.path(
    system.file(package = "joinme"),
    "stan",
    "joinme_mix_dynpred_threading.stan"
  )
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stan_file <- testthat::test_path(
      "..",
      "..",
      "inst",
      "stan",
      "joinme_mix_dynpred_threading.stan"
    )
  }
  required_data <- .stan_data_names(stan_file)
  expect_true(all(c(
    "use_dynamic_mixture",
    "dynamic_mix_probability_subject",
    "dynamic_mix_probability_marker",
    "dynamic_mix_location",
    "dynamic_mix_scale",
    "dynamic_mix_idx_subject",
    "dynamic_mix_idx_covariance",
    "dynamic_mix_idx_marker"
  ) %in% required_data))

  source <- paste(.read_stan_with_includes(stan_file), collapse = "\n")
  expect_match(source, "latent_progress_component_lpdf", fixed = TRUE)
  expect_match(
    source,
    "posterior_class_probability_new_subject",
    fixed = TRUE
  )
  expect_match(
    source,
    "posterior_class_probability_new_marker",
    fixed = TRUE
  )
})

test_that("ordinary dynamic prediction has no latent-class contract", {
  stan_file <- testthat::test_path(
    "..",
    "..",
    "inst",
    "stan",
    "joinme_dynpred_threading.stan"
  )
  required_data <- .stan_data_names(stan_file)
  source <- paste(.read_stan_with_includes(stan_file), collapse = "\n")

  expect_false(any(grepl("^dynamic_mix|^use_dynamic_mixture$", required_data)))
  expect_false(grepl("latent_progress_component_lpdf", source, fixed = TRUE))
  expect_false(
    grepl("posterior_class_probability_new_", source, fixed = TRUE)
  )
})

test_that("mixture arrays retain their Stan dimensions", {
  converted <- .coerce_rstan_mixture_data(list(
    mix_idx_subject = integer(0),
    mix_idx_marker = c(1L, 2L),
    mix_idx_covariance = integer(0),
    mix_probability_prior = c(1, 1, 1)
  ))

  expect_identical(dim(converted$mix_idx_subject), 0L)
  expect_identical(dim(converted$mix_idx_marker), 2L)
  expect_identical(dim(converted$mix_probability_prior), 3L)
})

test_that("prediction draw subsetting preserves mixture coordinate layouts", {
  draw_list <- list(
    beta_fixed = matrix(1:4, nrow = 2L),
    dynamic_mix_probability = matrix(c(0.4, 0.6, 0.7, 0.3), nrow = 2L),
    dynamic_mix_idx_subject = c(2L, 4L),
    dynamic_mix_idx_covariance = c(1L, 3L),
    dynamic_n_clusters = 2L
  )

  subset <- .subset_draws_for_prediction(draw_list, 1L)

  expect_identical(subset$dynamic_mix_idx_subject, c(2L, 4L))
  expect_identical(subset$dynamic_mix_idx_covariance, c(1L, 3L))
  expect_identical(subset$dynamic_n_clusters, 2L)
  expect_identical(nrow(subset$dynamic_mix_probability), 1L)
})

test_that("three-dimensional dynamic class draws retain paired uncertainty", {
  variable_names <- as.vector(sapply(1:2, function(fitted_draw) {
    as.vector(sapply(1:2, function(marker) {
      paste0(
        "posterior_class_probability_new_marker[",
        fitted_draw,
        ",",
        marker,
        ",",
        1:2,
        "]"
      )
    }))
  }))
  draw_matrix <- matrix(
    seq_len(2L * length(variable_names)),
    nrow = 2L,
    dimnames = list(NULL, variable_names)
  )

  extracted <- .extract_array3_from_stan(
    draw_matrix,
    "posterior_class_probability_new_marker",
    second_dimension = 2L,
    third_dimension = 2L,
    target_draws = 4L
  )

  expect_identical(dim(extracted), c(4L, 2L, 2L))
  expect_true(all(is.finite(extracted)))
  expect_equal(
    extracted[, 1L, 1L],
    unname(c(
      draw_matrix[1L, "posterior_class_probability_new_marker[1,1,1]"],
      draw_matrix[1L, "posterior_class_probability_new_marker[2,1,1]"],
      draw_matrix[2L, "posterior_class_probability_new_marker[1,1,1]"],
      draw_matrix[2L, "posterior_class_probability_new_marker[2,1,1]"]
    ))
  )
})

test_that("dynamic prediction exposes conditional class probabilities", {
  marker_probability <- array(
    c(
      0.8, 0.2, 0.7, 0.3,
      0.3, 0.7, 0.4, 0.6
    ),
    dim = c(2L, 2L, 2L)
  )
  prediction <- structure(
    list(
      draws = list(
        posterior_class = list(
          patient_1 = list(
            subject = matrix(
              c(0.9, 0.1, 0.8, 0.2),
              nrow = 2L,
              byrow = TRUE
            ),
            marker = marker_probability
          )
        )
      ),
      mixture = list(n_clusters = 2L),
      metadata = list(marker_levels = c("marker_a", "marker_b"))
    ),
    class = c("JoiNMeMixDynPred", "JoiNMeDynPred")
  )

  membership <- posterior_class(prediction)

  expect_s3_class(membership, "JoiNMePosteriorClass")
  expect_identical(
    unique(membership$subject$assigned_class),
    "class_1"
  )
  expect_setequal(
    unique(membership$marker$unit),
    c("marker_a", "marker_b")
  )
  expect_s3_class(
    plot.JoiNMeMixDynPred(prediction, type = "class_membership"),
    "ggplot"
  )
})

test_that("covariance class reconstruction returns a valid covariance matrix", {
  covariance <- .mixture_covariance_from_predictor(
    predictor = c(log(2), atanh(0.4), log(3)),
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )

  expect_equal(diag(covariance), c(4, 9), tolerance = 1e-10)
  expect_equal(
    covariance[1, 2] / sqrt(covariance[1, 1] * covariance[2, 2]),
    0.4,
    tolerance = 1e-10
  )
  expect_true(all(eigen(covariance, symmetric = TRUE)$values > 0))
})

test_that("marginal block draws retain component and ordinary coordinates", {
  centre <- array(0, dim = c(4L, 2L, 3L))
  component_draw <- array(
    seq_len(4L * 2L * 1L),
    dim = c(4L, 2L, 1L)
  )
  mixture <- list(
    dimensions = list(subject = 2L),
    starts = c(subject = 1L)
  )

  first <- .mixture_full_block_sample(
    block_centre = centre,
    sampled_latent = component_draw,
    mixture = mixture,
    level = "subject",
    seed = 47
  )
  second <- .mixture_full_block_sample(
    block_centre = centre,
    sampled_latent = component_draw,
    mixture = mixture,
    level = "subject",
    seed = 47
  )

  expect_equal(first, second)
  expect_equal(first[, , 2L], component_draw[, , 1L])
  expect_true(any(first[, , 1L] != 0))
  expect_true(any(first[, , 3L] != 0))
})

test_that("Cholesky reconstruction matches the covariance helper", {
  predictor <- c(log(1.5), atanh(-0.25), log(0.8))
  cholesky <- .mixture_cholesky_from_predictor(
    predictor,
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )
  covariance <- .mixture_covariance_from_predictor(
    predictor,
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )

  expect_equal(covariance, tcrossprod(cholesky))
  expect_equal(
    .assoc_corr_features_from_chol(cholesky),
    -0.25,
    tolerance = 1e-10
  )
})

test_that("class trajectories distinguish centres and response-scale margins", {
  fitted_model <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = list(
        P = 1L,
        R_id = 1L,
        R_mk = 1L,
        Q_idm = 1L,
        D = 1L,
        M_cov = 1L,
        K_cov = 0L,
        marker_levels = "marker_a",
        zidm_cols = "intercept",
        indep_idmarker_cov = 1L,
        allow_marker_crosscorr = 0L,
        vcov_diag_link = 1L,
        tmax = 1,
        idx_time_idm = integer(0),
        link_long = 1L,
        inv_link_n_ops = 0L,
        inv_link_n_const = 0L,
        inv_link_ops = matrix(integer(0), 1L, 0L),
        inv_link_const = matrix(numeric(0), 1L, 0L),
        Xcov = matrix(numeric(0), 1L, 0L)
      ),
      mixture = list(
        n_clusters = 2L,
        cluster_type = "subject",
        dimensions = list(
          subject = 1L,
          marker = integer(0),
          corr = integer(0),
          vcov = integer(0)
        ),
        starts = c(
          subject = 1L,
          marker = 0L,
          corr = 0L,
          vcov = 0L
        ),
        total_dimension = 1L,
        distribution = "normal"
      ),
      config = list(),
      dataLong = data.frame(
        id = 1L,
        marker = factor("marker_a"),
        time = 0,
        y = 0
      ),
      call = quote(joinme_mix(
        y ~ time,
        dataLong,
        id_var = "id",
        marker_var = "marker",
        time_var = "time"
      )),
      tmax = 1
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )

  testthat::local_mocked_bindings(
    .JoiNMefit_longitudinal_design_matrices = function(...) {
      list(
        fixed = matrix(c(1, 1), ncol = 1L),
        subject = matrix(c(1, 1), ncol = 1L),
        marker = matrix(c(1, 1), ncol = 1L),
        subject_marker = matrix(c(1, 1), ncol = 1L)
      )
    },
    .get_draws_matrix = function(fit, variables, draws = NULL, seed = 1) {
      if (all(grepl("^beta_scaled", variables))) {
        return(matrix(c(0, 0), ncol = 1L, dimnames = list(NULL, variables)))
      }
      if (all(grepl("^alpha_L", variables))) {
        return(matrix(c(0, 0), ncol = 1L, dimnames = list(NULL, variables)))
      }
      if (all(grepl("^lambda_L", variables))) {
        return(matrix(c(0, 0), ncol = 1L, dimnames = list(NULL, variables)))
      }
      stop("Unexpected posterior variable")
    },
    .mixture_matrix_draws = function(
      object,
      variable,
      number_rows,
      number_columns,
      draws,
      seed
    ) {
      output <- array(0, dim = c(2L, number_rows, number_columns))
      if (variable == "class_mean_subject") {
        output[, 1L, 1L] <- -1
        output[, 2L, 1L] <- 1
      }
      if (variable %in% c("L_u", "L_v")) {
        output[, 1L, 1L] <- 1
      }
      if (variable == "mix_location") {
        output[, 1L, 1L] <- -1
        output[, 2L, 1L] <- 1
      }
      if (variable == "mix_scale") {
        output[] <- 0.25
      }
      output
    },
    .apply_inverse_link_matrix = function(eta, ...) eta,
    .package = "joinme"
  )

  centre <- .mixture_class_trajectory(
    fitted_model,
    estimand = "mean_per_class",
    draws = 2L,
    seed = 5L,
    longitudinal_times = c(0, 1),
    longitudinal_points = 2L,
    marginal_samples = 2L,
    ci_levels = 0.95
  )
  marginal <- .mixture_class_trajectory(
    fitted_model,
    estimand = "marginal_per_class",
    draws = 2L,
    seed = 5L,
    longitudinal_times = c(0, 1),
    longitudinal_points = 2L,
    marginal_samples = 4L,
    ci_levels = 0.95,
    retain_association_samples = TRUE
  )

  expect_equal(dim(centre$epred[[1L]]), c(2L, 2L))
  expect_equal(dim(marginal$epred[[1L]]), c(2L, 2L))
  expect_true(all(centre$epred[[1L]] < centre$epred[[2L]]))
  expect_length(marginal$association_parts[[1L]], 4L)
  expect_true(all(c("mean", "marker", "total") %in%
    names(marginal$association_parts[[1L]][[1L]])))
  expect_false(isTRUE(all.equal(
    marginal$epred[[1L]],
    centre$epred[[1L]]
  )))
})

test_that("covariance association plots use the hazard contribution scale", {
  object <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      mixture = list(n_clusters = 2L, cluster_type = "vcov"),
      config = list(),
      stan_data = list(
        Q_idm = 1L,
        indep_idmarker_cov = 1L,
        M_cov = 1L,
        K_cov = 0L,
        Xcov = matrix(numeric(0), 1L, 0L),
        assoc_vcov = 1L,
        vcov_diag_link = 1L,
        idx_time_idm = integer(0),
        tmax = 1
      )
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )
  cached_association <- list(
    coeff_draws = list(vcov = matrix(c(0.5, 0.5), ncol = 1L)),
    term_map = data.frame(term = "vcov[1]", variable = "alpha_vcov[1]")
  )

  testthat::local_mocked_bindings(
    .mixture_matrix_draws = function(
      object,
      variable,
      number_rows,
      number_columns,
      draws,
      seed
    ) {
      output <- array(0, dim = c(2L, number_rows, number_columns))
      if (variable == "class_mean_covariance") {
        output[, 1L, 1L] <- -0.5
        output[, 2L, 1L] <- 0.5
      }
      output
    },
    .get_draws_matrix = function(fit, variables, draws = NULL, seed = 1) {
      value <- if (all(grepl("^lambda_L", variables))) 1 else 0
      matrix(
        value,
        nrow = 2L,
        ncol = length(variables),
        dimnames = list(NULL, variables)
      )
    },
    .get_association_plot_data = function(...) cached_association,
    .mixture_subset_association_data = function(
      association_data,
      ...
    ) association_data,
    .mixture_transform_draw_values = function(
      raw_values,
      ...
    ) raw_values,
    .association_transform_matrix = function(
      n_draws,
      ...
    ) matrix(0, nrow = n_draws, ncol = 1L),
    .package = "joinme"
  )

  association_plot_object <-
    .plot_mixture_covariance_association(
      object,
      association_term = "vcov[1]",
      draws = 2L,
      seed = 1L,
      ci_levels = 0.95,
      theme_fn = ggplot2::theme_bw
    )

  expect_s3_class(association_plot_object, "ggplot")
  expect_true(all(c("class", "x", "median") %in% names(
    association_plot_object$data
  )))
  expect_equal(
    length(unique(association_plot_object$data$class)),
    2L
  )
})

test_that("update preserves the mixture entry point and longitudinal-only mode", {
  data_for_update <- data.frame(
    id = 1L,
    marker = "marker_a",
    time = 0,
    y = 1
  )
  object <- structure(
    list(
      call = quote(joinme_mix(
        formulaLong = y ~ time,
        dataLong = data_for_update,
        n_clusters = 2,
        cluster_type = "subject",
        fit = FALSE
      )),
      formulaLong = y ~ time,
      formulaEvent = survival::Surv(
        .joinme_mix_follow_up,
        .joinme_mix_event
      ) ~ 1,
      formulaVCov = ~1,
      mixture = list(include_survival = FALSE),
      config = list()
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )
  captured <- NULL

  testthat::local_mocked_bindings(
    joinme_mix = function(...) {
      captured <<- list(...)
      "updated mixture"
    },
    .package = "joinme"
  )

  result <- update.JoiNMeFit(
    object,
    draws = 25L,
    .env = environment()
  )

  expect_identical(result, "updated mixture")
  expect_identical(captured$n_clusters, 2)
  expect_identical(captured$cluster_type, "subject")
  expect_identical(captured$draws, 25L)
  expect_null(captured$formulaEvent)
  expect_null(captured$dataEvent)
})
