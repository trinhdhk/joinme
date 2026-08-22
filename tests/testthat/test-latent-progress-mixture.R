test_that("ordinary models receive a parameter-free mixture representation", {
  mixture_data <- .build_mixture_standata(
    stan_data = list(),
    mixture = NULL
  )

  expect_identical(mixture_data$use_mixture, 0L)
  expect_identical(mixture_data$n_classes, 1L)
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
    marker_weight_sets_shared = 1L
  )
  mixture_data <- .build_mixture_standata(
    stan_data,
    mixture = list(
      n_classes = 3L,
      class_type = c("subject", "vcov"),
      class_dimensions = NULL,
      class_prior = list(baseline_prob = 1, slope = prior_normal(), family = prior_student_t(df = 6)),
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$n_classes, 3L)
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
      indep_idmarker_cov = 0L
    ),
    mixture = list(
      n_classes = 2L,
      class_type = "corr",
      class_dimensions = c(1L, 2L),
      class_prior = list(baseline_prob = 1, slope = prior_normal(), family = prior_normal()),
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$mix_idx_covariance, c(2L, 4L))
  expect_identical(mixture_data$mixture$class_type, "corr")
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
    marker_weight_sets_shared = 1L
  )
  mixture_data <- .build_mixture_standata(
    stan_data,
    mixture = list(
      n_classes = 2L,
      class_type = c("subject", "marker"),
      class_dimensions = NULL,
      class_prior = list(baseline_prob = c(2, 3), slope = prior_normal(), family = prior_normal()),
      include_survival = TRUE
    )
  )

  expect_identical(mixture_data$n_classes, 2L)
  expect_identical(mixture_data$K_mix, 4L)
  expect_equal(mixture_data$mix_probability_prior, c(2, 3))
  expect_identical(
    mixture_data$mixture$allocation_domains,
    list(subject = "subject", marker = "marker")
  )
})

test_that("class-regression priors are assembled in subject-then-marker order", {
  stan_data <- list(
    n_id = 3L,
    D = 2L,
    R_id = 2L,
    R_mk = 2L,
    Q_idm = 1L,
    indep_idmarker_cov = 1L
  )
  class_design <- list(
    subject = list(
      matrix = matrix(0, nrow = 3L, ncol = 2L),
      class_term_start = c(1L, 2L),
      class_term_count = c(1L, 1L)
    ),
    marker = list(
      matrix = matrix(0, nrow = 2L, ncol = 1L),
      class_term_start = c(1L, 1L),
      class_term_count = c(1L, 0L)
    )
  )
  specification <- list(
    n_classes = 2L,
    class_type = c("subject", "marker"),
    class_prior = list(
      baseline_prob = 1,
      slope = prior_laplace(
        mu = c(0, 0.5, -0.5),
        scale = c(1, 0.75, 0.5)
      ),
      family = prior_normal()
    ),
    class_design = class_design,
    include_survival = TRUE
  )

  mixture_data <- .build_mixture_standata(stan_data, specification)
  expect_identical(mixture_data$prior_class_regression_family, 3L)
  expect_equal(mixture_data$prior_class_regression_mu, c(0, 0.5, -0.5))
  expect_equal(mixture_data$prior_class_regression_scale, c(1, 0.75, 0.5))

  specification$class_prior$slope <- prior_normal(scale = c(1, 2))
  expect_error(
    .build_mixture_standata(stan_data, specification),
    "exactly one value per parameter"
  )
})

test_that("class_type accepts only the four public types", {
  expect_identical(
    .canonical_class_types(c("subject", "marker", "corr")),
    c("subject", "marker", "corr")
  )
  expect_error(.canonical_class_types("individual"), "subject")
  expect_error(.canonical_class_types("id_marker"), "subject")
  expect_error(.canonical_class_types("marker_weight"), "subject")
  expect_error(
    .canonical_class_types(c("corr", "vcov")),
    "cannot contain both"
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
    n_marker_weight_sets = 0L
  )
  specification <- list(
    n_classes = 2L,
    class_type = "subject",
    class_dimensions = 2L,
    class_prior = list(baseline_prob = 1, slope = prior_normal(), family = prior_student_t(df = 6)),
    include_survival = FALSE
  )

  expect_error(
    .build_mixture_standata(stan_data, specification),
    "Choose distinct integer indices between 1 and 1"
  )
  specification$class_type <- "marker"
  specification$class_dimensions <- NULL
  expect_error(
    .build_mixture_standata(stan_data, specification),
    "unavailable"
  )
  specification$class_type <- "subject"
  specification$n_classes <- 2.5
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
  scaffold <- .longitudinal_only_event_scaffold(
    longitudinal_data,
    id_variable = "id",
    time_variable = "time"
  )

  expect_equal(scaffold$id, c("a", "b"))
  expect_equal(scaffold$.joinme_follow_up, c(2, 2))
  expect_true(all(scaffold$.joinme_event == 0L))
})

test_that("formulaClass creates unit-level designs and stable probabilities", {
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
    n_classes = 3L
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
    .resolve_class_ordering(NULL, design$class_specific),
    "none"
  )
  expect_identical(
    .resolve_class_ordering("probability", design$class_specific),
    "probability"
  )
})

test_that("class_ordering targets only an intercept, baseline probabilities, or neither", {
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
    n_marker_weight_sets = 0L
  )
  specification <- list(
    n_classes = 3L,
    class_type = "subject",
    class_dimensions = c(1L, 2L),
    class_prior = list(baseline_prob = 1, slope = prior_normal(), family = prior_student_t(df = 6)),
    class_ordering = "intercept",
    include_survival = FALSE
  )

  intercept_ordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(intercept_ordered$mix_ordering, 1L)
  expect_identical(
    intercept_ordered$mix_ordered_location_coordinate,
    1L
  )

  specification$class_ordering <- "probability"
  probability_ordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(probability_ordered$mix_ordering, 2L)
  expect_identical(
    probability_ordered$mix_ordered_location_coordinate,
    0L
  )

  specification$class_ordering <- "none"
  unordered <- .build_mixture_standata(stan_data, specification)
  expect_identical(unordered$mix_ordering, 0L)
  expect_identical(unordered$mix_ordered_location_coordinate, 0L)
})

test_that("latent-class priors are declared through jm_prior", {
  priors <- jm_prior(
    class = list(
      baseline_prob = c(2, 3, 4),
      slope = prior_student_t(
        df = 4,
        mu = c(0, 0.5),
        scale = c(1.25, 0.75)
      ),
      family = prior_laplace()
    )
  )
  expect_s3_class(priors, "joinme_priors")
  expect_equal(priors$class$baseline_prob, c(2, 3, 4))
  expect_identical(priors$class$slope$family, "student_t")
  expect_equal(priors$class$slope$df, 4)
  expect_equal(priors$class$slope$mu, c(0, 0.5))
  expect_identical(priors$class$family$family, "laplace")
  expect_true("class" %in% names(formals(jm_prior)))
  expect_false("class_probability" %in% names(formals(jm_prior)))
  expect_false("class_regression" %in% names(formals(jm_prior)))
  expect_error(jm_prior(class = list(probability = 1)), "baseline_prob")
  expect_error(jm_prior(class = list(baseline_prob = 0)), "positive")
  expect_error(jm_prior(class = list(family = "normal")), "prior_normal")
  expect_error(
    jm_prior(class = list(family = prior_normal(mu = 0.25))),
    "standardised latent block"
  )
  expect_error(
    jm_prior(class = list(family = prior_laplace(scale = 1.5))),
    "standardised latent block"
  )
  expect_error(jm_prior(class = list(family = prior_horseshoe())), "does not accept")
  inherited <- jm_prior(slope = prior_laplace(scale = 0.75))
  expect_identical(inherited$class$slope$family, "laplace")
  expect_equal(inherited$class$slope$scale, 0.75)
  expect_false("class_probability_prior" %in% names(formals(joinme_mix)))
  expect_false("class_location_scale" %in% names(formals(joinme_mix)))
  expect_false("class_scale_rate" %in% names(formals(joinme_mix)))
  expect_false("class_separation" %in% names(formals(joinme_mix)))
})

test_that("joinme_mix prepares the inherited holder through the common entry path", {
  captured_arguments <- NULL
  fake_holder <- new.env(parent = emptyenv())
  fake_holder$stan_data <- list(mixture = list(
    n_classes = 2L,
    class_type = "subject",
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
    n_classes = 2,
    class_type = "subject",
    fit = FALSE
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(prepared$result_class, "mixture")
  expect_identical(captured_arguments$fit, FALSE)
  expect_identical(captured_arguments$mixture$n_classes, 2L)
  expect_identical(captured_arguments$mixture$class_type, "subject")
  expect_identical(captured_arguments$mixture$class_ordering, "intercept")
})

test_that("joinme_mix defaults class-specific formulae to no ordering", {
  captured_arguments <- NULL
  fake_holder <- new.env(parent = emptyenv())
  fake_holder$stan_data <- list(mixture = list(
    n_classes = 2L,
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
    n_classes = 2,
    class_type = "subject",
    formulaClass = list(~x, ~x),
    fit = FALSE
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(captured_arguments$mixture$class_ordering, "none")
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
    "n_classes",
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
    "prior_class_regression_family",
    "prior_class_regression_mu",
    "prior_class_regression_scale",
    "prior_class_regression_df",
    "prior_class_regression_global_df",
    "prior_class_regression_global_scale",
    "prior_class_regression_slab_df",
    "prior_class_regression_slab_scale"
  ) %in% required_data))

  source <- paste(.read_stan_with_includes(stan_file), collapse = "\n")
  master_source <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  expect_match(
    master_source,
    "#include include/submodels/latent_class/functions/component_density.stanfunctions",
    fixed = TRUE
  )
  expect_match(source, "include/submodels/latent_class/model/fit.stan", fixed = TRUE)
  expect_match(source, "ordered[n_classes] mix_location_ordered", fixed = TRUE)
  expect_match(source, "mix_location_unordered", fixed = TRUE)
  expect_match(
    source,
    "to_vector(mix_scale) ~ exponential(1)",
    fixed = TRUE
  )
  expect_match(source, "positive_ordered[n_classes - 1]", fixed = TRUE)
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
  master_source <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")

  expect_false(any(grepl("^mix_|^n_classes$|^K_mix$|^use_mixture$", required_data)))
  expect_false(grepl("include/submodels/latent_class/functions", master_source, fixed = TRUE))
  expect_false(grepl("latent_progress_component_lpdf", source, fixed = TRUE))
  expect_false(grepl("posterior_class_probability_", source, fixed = TRUE))
})

test_that("mixture masters import only the function definitions they use", {
  fit_master <- paste(
    readLines(testthat::test_path("..", "..", "inst", "stan", "joinme_mix_fit_threading.stan"), warn = FALSE),
    collapse = "\n"
  )
  prediction_master <- paste(
    readLines(testthat::test_path("..", "..", "inst", "stan", "joinme_mix_dynpred_threading.stan"), warn = FALSE),
    collapse = "\n"
  )

  component_include <- "include/submodels/latent_class/functions/component_density.stanfunctions"
  probability_include <- "include/submodels/latent_class/functions/class_probability.stanfunctions"
  expect_true(grepl(component_include, fit_master, fixed = TRUE))
  expect_true(grepl(probability_include, fit_master, fixed = TRUE))
  expect_true(grepl(component_include, prediction_master, fixed = TRUE))
  expect_false(grepl(probability_include, prediction_master, fixed = TRUE))
})

test_that("the mixture prior is evaluated once outside the threaded likelihood", {
  stan_file <- testthat::test_path(
    "..",
    "..",
    "inst",
    "stan",
    "joinme_mix_fit_threading.stan"
  )
  main_source <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  mixture_prior_position <- regexpr(
    "#include include/submodels/latent_class/model/fit.stan",
    main_source,
    fixed = TRUE
  )[[1L]]
  likelihood_position <- regexpr(
    "#include include/etc/model/fit_threaded_likelihood.stan",
    main_source,
    fixed = TRUE
  )[[1L]]

  expect_gt(mixture_prior_position, 0L)
  expect_gt(likelihood_position, mixture_prior_position)
  expect_equal(
    lengths(gregexpr(
      "#include include/submodels/latent_class/model/fit.stan",
      main_source,
      fixed = TRUE
    )),
    1L
  )

  partial_source <- paste(readLines(
    testthat::test_path(
      "..",
      "..",
      "inst",
      "stan",
      "include",
      "etc",
      "functions",
      "joinme_fit_partial.stanfunctions"
    ),
    warn = FALSE
  ), collapse = "\n")
  expect_false(grepl(
    "latent_progress_component_lpdf",
    partial_source,
    fixed = TRUE
  ))
  expect_false(grepl("mix_scale", partial_source, fixed = TRUE))
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
    dynamic_n_classes = 2L
  )

  subset <- .subset_draws_for_prediction(draw_list, 1L)

  expect_identical(subset$dynamic_mix_idx_subject, c(2L, 4L))
  expect_identical(subset$dynamic_mix_idx_covariance, c(1L, 3L))
  expect_identical(subset$dynamic_n_classes, 2L)
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
      mixture = list(n_classes = 2L),
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
  expect_named(
    membership$subject,
    c(
      "unit",
      "class",
      "Estimate",
      "Est.Error",
      "Q2.5",
      "Q97.5",
      "Rhat",
      "ess_bulk",
      "ess_tail",
      "assigned_class"
    )
  )
  expect_s3_class(
    plot(membership, domain = "subject"),
    "ggplot"
  )
  printed_membership <- paste(
    capture.output(print(membership, max_rows = 2L)),
    collapse = "\n"
  )
  expect_match(
    printed_membership,
    "Maximum-probability allocation counts",
    fixed = TRUE
  )
  expect_match(
    printed_membership,
    "Estimate",
    fixed = TRUE
  )
  expect_s3_class(
    plot.JoiNMeMixDynPred(prediction, type = "class_membership"),
    "ggplot"
  )
})

test_that("class association relevance follows the fitted random-effect block", {
  subject_fit <- structure(
    list(
      mixture = list(class_type = "subject"),
      config = list()
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )
  covariance_fit <- structure(
    list(
      mixture = list(class_type = "vcov"),
      config = list()
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )

  expect_true(.mixture_class_changes_association(
    subject_fit,
    "cv_mean",
    "mean_per_class"
  ))
  expect_true(.mixture_class_changes_association(
    subject_fit,
    "cs_total",
    "mean_per_class"
  ))
  expect_false(.mixture_class_changes_association(
    subject_fit,
    "cv_marker",
    "mean_per_class"
  ))
  expect_true(.mixture_class_changes_association(
    covariance_fit,
    "vcov[1]",
    "mean_per_class"
  ))
  expect_false(.mixture_class_changes_association(
    covariance_fit,
    "cv_mean",
    "mean_per_class"
  ))
})

test_that("posterior marker weights retain draw-wise trajectory pairing", {
  trajectory <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3L,
    byrow = TRUE
  )
  weight <- posterior::as_draws_matrix(matrix(
    c(2, 3, 4),
    ncol = 1L,
    dimnames = list(NULL, "weight")
  ))

  weighted <- .mixture_weight_draw_trajectory(
    trajectory,
    weight[, 1L]
  )

  expect_equal(
    weighted,
    trajectory * c(2, 3, 4)
  )
})

test_that("covariance class reconstruction returns a valid covariance matrix", {
  covariance <- .mixture_covariance_from_predictor(
    predictor = c(log(2), atanh(0.4), log(3)),
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )

  expect_equal(diag(covariance), c(4, 9 * (1 + 0.4^2)), tolerance = 1e-10)
  expect_equal(
    covariance[1, 2] / covariance[2, 2]^0.5 / covariance[1, 1]^0.5,
    0.4 / sqrt(1 + 0.4^2),
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
  expect_equal(cholesky[2, 2] / exp(log(0.8)), 1, tolerance = 1e-12)
  expect_equal(
    .assoc_corr_features_from_chol(cholesky),
    -0.25 / sqrt(1 + 0.25^2),
    tolerance = 1e-10
  )
})

test_that("independent SD and K regression lets dependence change marginal variance", {
  predictor_zero_corr <- c(log(2), atanh(0.0), log(3))
  predictor_high_corr <- c(log(2), atanh(0.8), log(3))

  covariance_zero_corr <- .mixture_covariance_from_predictor(
    predictor = predictor_zero_corr,
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )
  covariance_high_corr <- .mixture_covariance_from_predictor(
    predictor = predictor_high_corr,
    q_dimension = 2L,
    diagonal_only = FALSE,
    diagonal_link = 1L
  )

  # In the independent SD+K parameterisation, off-diagonal K entries contribute
  # to marginal variance through K K^T.
  expect_equal(diag(covariance_zero_corr), c(4, 9), tolerance = 1e-10)
  expect_equal(diag(covariance_high_corr), c(4, 9 * (1 + 0.8^2)), tolerance = 1e-10)
  expect_gt(abs(covariance_high_corr[1, 2]), abs(covariance_zero_corr[1, 2]))
})

test_that("mixture summaries reduce membership to active class totals", {
  subject_membership <- data.frame(
    unit = rep(c("subject_a", "subject_b"), each = 3L),
    class = rep(paste0("class_", 1:3), times = 2L),
    probability = c(0.8, 0.1, 0.1, 0.2, 0.7, 0.1),
    lower = 0,
    upper = 1,
    assigned_class = rep(c("class_1", "class_2"), each = 3L),
    stringsAsFactors = FALSE
  )
  subject_draws <- matrix(
    c(
      0.8, 0.2, 0.1, 0.7, 0.1, 0.1,
      0.6, 0.4, 0.2, 0.5, 0.2, 0.1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      draw = 1:2,
      variable = c(
        "probability[1,1]",
        "probability[2,1]",
        "probability[1,2]",
        "probability[2,2]",
        "probability[1,3]",
        "probability[2,3]"
      )
    )
  ) # two posterior draws ordered by class and then allocation unit

  overview <- .mixture_membership_overview(
    class_membership = list(subject = subject_membership),
    class_membership_draws = list(subject = subject_draws)
  )

  expect_named(overview, "subject")
  expect_equal(nrow(overview$subject), 3L)
  expect_named(
    overview$subject,
    c(
      "class",
      "Estimate",
      "Est.Error",
      "Q2.5",
      "Q97.5",
      "assigned_units"
    )
  )
  expect_equal(
    overview$subject$Estimate,
    c(1, 0.75, 0.25)
  )
  expect_equal(
    overview$subject$assigned_units,
    c(1L, 1L, 0L)
  )
  expect_false("marker" %in% names(overview))
})

test_that("mixture class-regression summary skips zero fitted dimensions", {
  object <- structure(
    list(
      mixture = list(
        class_design = list(
          subject = list(
            # Legacy serialized objects may retain this placeholder label even
            # when no subject-domain class-regression coefficient was fitted.
            columns = "class_:"
          )
        )
      ),
      config = list(),
      stan_data = list(P_class_subject = 0L),
      fit = NULL
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )

  expect_null(.mixture_class_regression_summary(
    object = object,
    domain = "subject",
    number_classes = 2L,
    draws = 20,
    seed = 1,
    digits = 3
  ))
})

test_that("mixture class-regression summary skips absent draw variables", {
  object <- structure(
    list(
      mixture = list(
        class_design = list(
          subject = list(
            coefficient_covariate = "treatment",
            coefficient_class = 1L
          )
        )
      ),
      config = list(),
      stan_data = list(P_class_subject = 1L),
      fit = NULL
    ),
    class = c("JoiNMeMixFit", "JoiNMeFit")
  )

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, ...) {
      posterior::as_draws_array(array(
        0,
        dim = c(2L, 1L, 1L),
        dimnames = list(
          iteration = c("1", "2"),
          chain = "1",
          variable = "other_parameter[1]"
        )
      ))
    },
    .package = "joinme"
  )

  expect_null(.mixture_class_regression_summary(
    object = object,
    domain = "subject",
    number_classes = 2L,
    draws = 20,
    seed = 1,
    digits = 3
  ))
})

test_that("repeated scientific labels retain one distinct posterior row per variable", {
  draw_array <- array(
    NA_real_,
    dim = c(20L, 2L, 4L),
    dimnames = list(
      iteration = seq_len(20L),
      chain = seq_len(2L),
      variable = paste0("coefficient[", seq_len(4L), "]")
    )
  ) # four distinct source variables whose display labels repeat by class
  for (variable_index in seq_len(4L)) {
    draw_array[, , variable_index] <- variable_index +
      matrix(seq_len(40L) / 1000, nrow = 20L, ncol = 2L)
  }

  summary_table <- .assoc_summary_from_draw_array(
    draw_array,
    term_labels = rep(c("treatment", "age"), times = 2L)
  )
  mixture_table <- .mixture_summary_columns(
    summary_table,
    identifier_columns = "term"
  )

  expect_equal(nrow(mixture_table), 4L)
  expect_equal(mixture_table$term, c("treatment", "age", "treatment", "age"))
  expect_equal(length(unique(mixture_table$Estimate)), 4L)
  expect_named(
    mixture_table,
    c(
      "term",
      "Estimate",
      "Est.Error",
      "Q2.5",
      "Q97.5",
      "Rhat",
      "ess_bulk",
      "ess_tail"
    )
  )
})

test_that("mixture headline diagnostics retain sampler-wide extremes", {
  displayed_parameters <- data.frame(
    term = c("class_1", "class_2"),
    Rhat = c(1.02, 1.08),
    ess_bulk = c(80, 140),
    ess_tail = c(120, 70)
  ) # two displayed class parameters used for threshold counts
  sampler_diagnostics <- list(
    draws = 1000,
    divergences = 2,
    treedepth_hits = 3,
    ebfmi_min = 0.12,
    max_rhat = 1.15,
    min_ess_bulk = 40,
    min_ess_tail = 35
  ) # sampler-wide diagnostics, including variables absent from the table

  diagnostic_table <- .summary_diagnostics_table(
    sampler_diagnostics = sampler_diagnostics,
    reported_tables = list(class_location = displayed_parameters)
  )
  diagnostic_value <- stats::setNames(
    diagnostic_table$value,
    diagnostic_table$metric
  ) # named vector for direct assertions on the common diagnostic schema

  expect_equal(diagnostic_value[["max_rhat"]], 1.15)
  expect_equal(diagnostic_value[["min_ess_bulk"]], 40)
  expect_equal(diagnostic_value[["min_ess_tail"]], 35)
  expect_equal(diagnostic_value[["n_terms_total"]], 2)
  expect_equal(diagnostic_value[["n_terms_bad_rhat"]], 2)
  expect_equal(diagnostic_value[["n_terms_low_ess_bulk"]], 1)
  expect_equal(diagnostic_value[["n_terms_low_ess_tail"]], 1)
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
        K_cov_sd = 0L,
        K_cov_corr = 0L,
        marker_levels = "marker_a",
        zidm_cols = "intercept",
        indep_idmarker_cov = 1L,
        allow_marker_crosscorr = 0L,
        vcov_diag_link = 1L,
        tmax = 1,
        link_long = 1L,
        inv_link_n_ops = 0L,
        inv_link_n_const = 0L,
        inv_link_ops = matrix(integer(0), 1L, 0L),
        inv_link_const = matrix(numeric(0), 1L, 0L),
        Xcov_sd = matrix(numeric(0), 1L, 0L),
        Xcov_corr = matrix(numeric(0), 1L, 0L)
      ),
      mixture = list(
        n_classes = 2L,
        class_type = "subject",
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
    ..longitudinal_design_matrices = function(...) {
      list(
        fixed = matrix(c(1, 1), ncol = 1L),
        subject = matrix(c(1, 1), ncol = 1L),
        marker = matrix(c(1, 1), ncol = 1L),
        subject_marker = matrix(c(1, 1), ncol = 1L)
      )
    },
    .get_draws_matrix = function(fit, variables, draws = NULL, seed = 1) {
      if (all(grepl("^beta\\[", variables))) {
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
      mixture = list(n_classes = 2L, class_type = "vcov"),
      config = list(),
      stan_data = list(
        Q_idm = 1L,
        indep_idmarker_cov = 1L,
        M_cov = 1L,
        K_cov_sd = 0L,
        K_cov_corr = 0L,
        Xcov_sd = matrix(numeric(0), 1L, 0L),
        Xcov_corr = matrix(numeric(0), 1L, 0L),
        assoc_vcov = 1L,
        vcov_diag_link = 1L,
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
        n_classes = 2,
        class_type = "subject",
        fit = FALSE
      )),
      formulaLong = y ~ time,
      formulaEvent = survival::Surv(
        .joinme_follow_up,
        .joinme_event
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
  expect_identical(captured$n_classes, 2)
  expect_identical(captured$class_type, "subject")
  expect_identical(captured$draws, 25L)
  expect_null(captured$formulaEvent)
  expect_null(captured$dataEvent)
})
