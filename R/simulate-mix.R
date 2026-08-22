#' Simulate data from a latent-progress JoiNMe model
#'
#' @description
#' `simulate_joinme_mix()` is the generative counterpart of [joinme_mix()].
#' It follows the ordinary [simulate_joinme()] workflow, but replaces selected
#' standardised random-effect coordinates by a finite location--scale mixture.
#' The selected coordinates, class-regression formulae, ordering convention
#' and shared-allocation rules have the same meaning as in [joinme_mix()].
#'
#' The function deliberately delegates observation-family sampling,
#' distributional regression, marker-weight construction, event-time
#' generation, censoring, delayed entry and association transformations to
#' [simulate_joinme()]. Only the distribution of the selected standardised
#' latent coordinates is changed. This provides a direct simulation-to-fit
#' workflow without maintaining a second longitudinal or survival simulator.
#' Population parameters follow the same rule as [simulate_joinme()]: numeric
#' declarations in `truth` are fixed exactly and `prior_*()` declarations are
#' drawn once. The realised values are held fixed for the simulated data set
#' and stored in `truth`. This includes shared or association-term-specific
#' marker-weight locations. Class probabilities, regression coefficients,
#' locations and scales are recorded under `truth$mixture`.
#'
#' @details
#' Let \eqn{C_j\in\{1,\ldots,G\}} denote the latent class for allocation unit
#' \eqn{j}. Conditional on class \eqn{g}, every selected coordinate \eqn{r} is drawn
#' independently as
#'
#' \deqn{
#' z_{jr}\mid C_j=g \sim
#' D\{\mu_{gr},\sigma_{gr}\},
#' }
#'
#' where \eqn{D} is declared by
#' `jm_truth(class = list(family = prior_student_t(...)))`,
#' [prior_normal()], or [prior_laplace()]. The same declaration is retained as
#' the recovery prior used by [joinme_mix()].
#' Unselected subject, marker and covariance-regression coordinates are
#' standard Normal.
#'
#' Marker weights are not latent-class coordinates. When they are estimated,
#' their common location is supplied by `truth$marker_weights$intercept`,
#' and direct marker departures are drawn from the centred unit-scale
#' `truth$marker_weights$family`. No further scale multiplies those
#' departures because the association slope already scales the weighted
#' marker feature. The family name `"student_t"` draws one value shared by all
#' weight sets as `2 + Gamma(2, 0.1)`. The realised value is retained in `truth`. The departure
#' coordinate retains location zero and ordinary scale one. A constant family
#' bypasses both quantities, so the effective weights equal the offsets declared
#' in `truth$marker_weights` exactly.
#'
#' Subject and covariance-regression coordinates share one subject-level
#' allocation. Combining compatible types consequently retains exactly `n_classes`
#' components rather than constructing a Cartesian product.
#'
#' `formulaClass` is evaluated using the same allocation-unit design builder
#' as [joinme_mix()]. A shared formula is applied to every active allocation
#' domain. A named list may supply separate `subject` and `marker`
#' formulae. A list of exactly `n_classes` formulae supplies class-specific
#' predictors. The formula intercept is omitted because `probability` in
#' `class_parameters` already specifies the zero-covariate probabilities.
#'
#' `class_parameters` is a named list with the following entries:
#'
#' - `probability`: a positive vector of length `n_classes`, normalised to sum
#'   to one. If omitted, positive values are drawn once and normalised; under
#'   probability ordering the realised probabilities are sorted increasingly.
#' - `location`: an `n_classes` by \(K\) matrix in the packed coordinate
#'   order, or a named list containing matrices for selected class
#'   types. Omitted locations are drawn once around separated class centres;
#'   the ordered coordinate is sorted when intercept ordering is requested.
#' - `scale`: a positive scalar, vector, matrix, or named level list matching
#'   `location`. Omitted component scales are drawn once from `Uniform(0.45,
#'   0.85)`.
#'
#' The returned `truth$mixture` record contains the realised class
#' probabilities, allocations, component parameters, coordinate layout and
#' standardised latent draws. Its `stan_data_contract` entry records the exact
#' mixture dimensions, source indices, packed starts, ordering code and
#' class-design segments expected from a corresponding `joinme_mix(...,
#' fit = FALSE)` call. Class labels are kept out of `dataLong` and `dataEvent`
#' so that an analysis cannot accidentally use the simulated truth as an
#' observed predictor.
#'
#' A longitudinal-only simulation is requested with `formulaEvent = NULL`,
#' matching [joinme_mix()]. The common engine then uses an internal
#' zero-hazard scaffold solely to construct covariates and observation rows;
#' the returned `dataEvent` and survival formula are `NULL`. Explicit
#' association terms are not permitted in this mode.
#'
#' The `truth` argument is created with [jm_truth()] and follows the scientific
#' component hierarchy used for fitting. It describes data-generating
#' quantities rather than fitting distributions. In particular,
#' class truths are declared together as
#' `class = list(baseline_prob = ..., slope = ..., family = ...)`:
#' `baseline_prob` governs the Dirichlet prior, `slope` governs coefficients
#' introduced by `formulaClass`, and `family` governs the centred unit-scale
#' distribution within every latent class. A numeric `class$slope` fixes the generating coefficient
#' vector; a `prior_*()` declaration draws it once at the beginning of the
#' simulation. Class probabilities, locations and scales are supplied through
#' `class_parameters` when fixed generating values are required.
#'
#' @inheritParams simulate_joinme
#' @param n_classes Number of latent classes \(G\); must be an integer of at
#'   least two. This has the same meaning as in [joinme_mix()].
#' @param formulaClass One-sided shared, domain-specific, or class-specific
#'   class-membership formula, following [joinme_mix()].
#' @param class_type Character vector selecting only `"subject"`, `"marker"`,
#'   `"corr"` and/or `"vcov"`, with the same meaning as in [joinme_mix()].
#' @param class_dimensions Optional coordinate selection with the same
#'   syntax as [joinme_mix()].
#' @param class_ordering Class-label identification rule: `"intercept"`,
#'   `"probability"` or `"none"`. The contextual default is identical to
#'   [joinme_mix()].
#' @param class_parameters Generative class probabilities, component locations
#'   and component scales. See Details.
#'
#' @return A list with the same top-level components as [simulate_joinme()].
#'   The `truth` entry additionally contains a `mixture` record suitable for
#'   recovery studies and class-allocation assessment. It also records the
#'   common marker-weight location as `marker_weight_mean`: one value named
#'   `shared` when weight sets are shared, or one named value per active
#'   association term otherwise. `marker_weight_mean_by_term` provides the
#'   corresponding association-term representation; entries for inactive
#'   weighted association forms are zero.
#'   `truth$recovery$entry_point` is `"joinme_mix"`, and its `arguments`
#'   preserve the formula-class, class-type, coordinate and ordering syntax
#'   alongside every ordinary fitting argument.
#'
#' @examples
#' \dontrun{
#' simulated <- simulate_joinme_mix(
#'   formulaLong = y ~ 1 + time + x1 +
#'     (1 + time | id) +
#'     (0 + x1 + (1 + time | id) | marker),
#'   formulaEvent = survival::Surv(time, event) ~ x1 + x2,
#'   n_id = 80,
#'   families = c("gaussian", "student_t", "gaussian"),
#'   n_classes = 3,
#'   formulaClass = ~ x1,
#'   class_type = c("subject", "vcov"),
#'   truth = jm_truth(class = list(slope = c(
#'     "class_1:x1" = -0.6,
#'     "class_2:x1" = 0.4
#'   ))),
#'   class_parameters = list(
#'     probability = c(0.25, 0.45, 0.30),
#'     location = matrix(
#'       c(-1.5, -0.8, 0, 0, 1.5, 0.8),
#'       nrow = 3,
#'       byrow = TRUE
#'     ),
#'     scale = 0.7
#'   ),
#'   seed = 701
#' )
#'
#' fitted <- joinme_mix(
#'   formulaLong = simulated$truth$formulaLong,
#'   dataLong = simulated$dataLong,
#'   formulaEvent = simulated$truth$formulaEvent,
#'   dataEvent = simulated$dataEvent,
#'   formulaVCov = simulated$truth$formulaVCov,
#'   families = simulated$marker_info$families,
#'   n_classes = 3,
#'   formulaClass = ~ x1,
#'   class_type = c("subject", "vcov")
#' )
#' }
#' @export
simulate_joinme_mix <- function(
  formulaLong = y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  formulaVCov = ~1,
  formulaDist = NULL,
  formulaAssoc = NULL,
  transforms = NULL,
  truth = joinme_truth(),
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  seed = .Random.seed[[1]],
  covariate_formulas = list(
    x1 ~ rnorm(n_id),
    x2 ~ rnorm(n_id)
  ),
  assoc = c("cv_total"),
  vcov_diag_link = "softplus",
  family_params = list(
    gaussian = list(sigma = 1.0),
    student_t = list(sigma = 1.5, nu = 4),
    bernoulli = list(),
    binomial = list(trials = 10),
    poisson = list(),
    negbin2 = list(phi = 2),
    skew_normal = list(sigma = 1.0, alpha = 0),
    double_exponential = list(sigma = 1.0),
    skew_double_exponential = list(sigma = 1.0, tau = 0.5),
    beta = list(kappa = 10),
    cumulative_logit = list(cutpoints = c(-1, 1))
  ),
  time_cens = 8.0,
  eps_cs = 1e-3,
  integration_control = list(
    rel.tol = 1e-6,
    subdivisions = 2000L,
    stop.on.error = TRUE
  ),
  root_control = list(
    t_init = 1.0,
    t_max = 50.0,
    expand = 1.7,
    max_expand = 60L
  ),
  quadrature_nodes = 15L,
  n_workers = 1L,
  use_mirai = TRUE,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  y_var = "y",
  event_time_var = "time",
  event_var = "event",
  n_classes = 2L,
  formulaClass = ~1,
  class_type = "subject",
  class_dimensions = NULL,
  class_ordering = NULL,
  class_parameters = list()
) {
  longitudinal_only <- is.null(
    formulaEvent
  ) # whether the public simulation omits the survival process
  supplied_call <- match.call(
    expand.dots = TRUE
  ) # original call retained to distinguish defaults from explicit requests
  if (!inherits(truth, "joinme_truth")) {
    cli::cli_abort(c(
      x = "{.arg truth} must be created with {.fn jm_truth}.",
      i = "Fixed generating values and their between-simulation distributions belong in that declaration."
    ))
  }
  selected_class_types <- .canonical_class_types(
    class_type
  ) # checked public class types before ordinary simulation work begins
  simulation_truth_request <- unclass(truth) # fixed-or-random population declarations used once for this mixture simulation
  checked_prior_specification <- attr(truth, "fitting_priors", exact = TRUE) # probability distributions retained for the optional recovery fit
  if (length(selected_class_types) == 0L) {
    cli::cli_abort(
      "{.arg class_type} must select at least one class type."
    )
  }
  if (
    longitudinal_only &&
      (
        "formulaAssoc" %in% names(supplied_call) &&
          !is.null(formulaAssoc)
      )
  ) {
    cli::cli_abort(
      "{.arg formulaAssoc} requires a survival process in {.fn simulate_joinme_mix}."
    )
  }
  if (
    longitudinal_only &&
      "assoc" %in% names(supplied_call) &&
      length(assoc) > 0L
  ) {
    cli::cli_abort(
      "{.arg assoc} requires a survival process in {.fn simulate_joinme_mix}."
    )
  }

  # Retain the fitting entry point's mixture syntax in one small,
  # unevaluated specification. The ordinary simulator will validate the
  # formula-derived dimensions after it has built the same random-effect
  # design matrices used for data generation.
  mixture_specification <- list(
    n_classes = n_classes, # requested number of latent components
    formulaClass = formulaClass, # covariate model for component membership
    class_type = selected_class_types, # random-effect blocks receiving mixture distributions
    class_dimensions = class_dimensions, # selected coordinates by block
    class_ordering = class_ordering, # label-identification convention
    class_parameters = class_parameters, # fixed generative mixture probabilities, locations and scales
    class_slope_declaration = simulation_truth_request$class$slope %||%
      checked_prior_specification$class$slope, # fixed vector or generating distribution for formulaClass coefficients
    class_prior = checked_prior_specification$class # recovery priors for baseline probabilities and class regression
  )

  # Reconstruct the ordinary simulation call from the arguments that
  # were actually supplied. Missing arguments remain missing, so their defaults
  # continue to be owned by `simulate_joinme()` rather than duplicated at
  # runtime. Removing only the six mixture arguments establishes a one-to-one
  # ordinary simulation path.
  simulation_arguments <- as.list(
    supplied_call
  )[-1L] # supplied ordinary and mixture arguments as language objects
  mixture_argument_names <- c(
    "n_classes",
    "formulaClass",
    "class_type",
    "class_dimensions",
    "class_ordering",
    "class_parameters"
  ) # arguments consumed by the mixture layer rather than the ordinary engine
  simulation_arguments[mixture_argument_names] <- NULL
  simulation_arguments <- lapply(
    simulation_arguments,
    eval,
    envir = parent.frame()
  ) # concrete values forwarded to the established ordinary simulator
  simulation_arguments$.mixture_specification <-
    mixture_specification # private, already checked mixture request

  # Execute the common generator. The private specification changes
  # only standardised random-effect draws; all downstream trajectory, response
  # and event calculations remain those of `simulate_joinme()`.
  simulation <- do.call(
    simulate_joinme,
    simulation_arguments,
    envir = parent.frame()
  )
  if (longitudinal_only) {
    simulation$dataEvent <- NULL
    simulation$truth$formulaEvent <- NULL
    simulation$truth$assoc <- character(0)
    simulation$truth$recovery$arguments$formulaEvent <- NULL
    simulation$truth$recovery$arguments$dataEvent <- NULL
    simulation$truth$recovery$arguments$assoc <- character(0)
  }
  simulation
}

#' Resolve packed mixture locations or scales
#'
#' @param supplied Scalar, vector, matrix, named level list, or `NULL`.
#' @param default Complete default matrix.
#' @param layout Checked coordinate layout from `.build_mixture_standata()`.
#' @param argument_name User-facing parameter name.
#' @param require_positive Whether every supplied value must be strictly
#'   positive.
#'
#' @return An `n_classes` by \(K\) numeric matrix.
#' @keywords internal
#' @noRd
.sim_mixture_component_matrix <- function(
  supplied,
  default,
  layout,
  argument_name,
  require_positive = FALSE
) {
  n_classes <- nrow(default) # number of component rows
  total_dimension <- ncol(default) # packed selected-coordinate count
  output <- default # complete matrix, overwritten only where requested

  normalise_block <- function(value, rows, columns, label) {
    if (is.matrix(value)) {
      value <- unname(value)
      if (!identical(dim(value), c(rows, columns))) {
        cli::cli_abort(c(
          x = "{.arg {argument_name}} for {label} has incompatible dimensions.",
          i = "Expected {rows} by {columns}; received {nrow(value)} by {ncol(value)}."
        ))
      }
      return(matrix(as.numeric(value), nrow = rows, ncol = columns))
    }
    value <- as.numeric(value)
    if (length(value) == 1L) {
      return(matrix(value, nrow = rows, ncol = columns))
    }
    if (length(value) == columns) {
      return(matrix(rep(value, each = rows), nrow = rows, ncol = columns))
    }
    if (columns == 1L && length(value) == rows) {
      return(matrix(value, nrow = rows, ncol = columns))
    }
    if (length(value) == rows * columns) {
      return(matrix(value, nrow = rows, ncol = columns, byrow = TRUE))
    }
    cli::cli_abort(c(
      x = "{.arg {argument_name}} for {label} has incompatible length.",
      i = "Supply a scalar, {columns} coordinate values, or a {rows} by {columns} matrix."
    ))
  }

  if (is.null(supplied)) {
    return(output)
  }
  if (is.list(supplied) && !is.matrix(supplied)) {
    if (is.null(names(supplied)) || any(!nzchar(names(supplied)))) {
      cli::cli_abort(
        "{.arg {argument_name}} must use class-type names when supplied as a list."
      )
    }
    canonical_names <- .canonical_class_types(names(supplied))
    if (anyDuplicated(canonical_names)) {
      cli::cli_abort(
        "{.arg {argument_name}} supplies the same class type more than once."
      )
    }
    names(supplied) <- canonical_names
    unknown_levels <- setdiff(
      canonical_names,
      layout$mixture$class_type
    )
    if (length(unknown_levels) > 0L) {
      cli::cli_abort(
        "{.arg {argument_name}} includes unselected type{?s}: {paste(unknown_levels, collapse = ', ')}."
      )
    }
    for (level in canonical_names) {
      block_dimension <- length(layout$mixture$dimensions[[level]])
      block_start <- layout$mixture$starts[[level]]
      block_columns <- seq.int(block_start, length.out = block_dimension)
      output[, block_columns] <- normalise_block(
        supplied[[level]],
        rows = n_classes,
        columns = block_dimension,
        label = level
      )
    }
  } else {
    output <- normalise_block(
      supplied,
      rows = n_classes,
      columns = total_dimension,
      label = "the packed mixture"
    )
  }

  if (any(!is.finite(output))) {
    cli::cli_abort("{.arg {argument_name}} must contain only finite values.")
  }
  if (isTRUE(require_positive) && any(output <= 0)) {
    cli::cli_abort("{.arg {argument_name}} must contain only positive values.")
  }
  output
}

#' Prepare a generative latent-progress mixture
#'
#' @param specification Private request assembled by `simulate_joinme_mix()`.
#' @param data_event One-row-per-subject simulation covariate data.
#' @param marker_prototype Marker-aligned prototype longitudinal data.
#' @param dimensions Named dimensions of the three latent blocks.
#' @param labels Named random-effect column labels.
#' @param covariance_is_diagonal Whether the marker-by-subject covariance is
#'   diagonal.
#' @param id_variable Subject identifier name.
#' @param marker_variable Marker identifier name.
#'
#' @return A fully checked layout with probabilities, allocations and
#'   component parameters.
#' @keywords internal
#' @noRd
.sim_prepare_mixture <- function(
  specification,
  data_event,
  marker_prototype,
  dimensions,
  labels,
  covariance_is_diagonal,
  id_variable,
  marker_variable
) {
  if (is.null(specification)) {
    return(NULL)
  }

  # Step 1: validate the public class count and selected class types
  # through the same helpers used by the fitting entry point.
  n_classes <- specification$n_classes # requested component count
  if (
    !is.numeric(n_classes) ||
      length(n_classes) != 1L ||
      !is.finite(n_classes) ||
      n_classes < 2 ||
      n_classes != floor(n_classes)
  ) {
    cli::cli_abort("{.arg n_classes} must be a single integer of at least two.")
  }
  n_classes <- as.integer(n_classes) # checked component count
  selected_types <- .canonical_class_types(
    specification$class_type
  ) # canonical selected random-effect blocks
  if (length(selected_types) == 0L) {
    cli::cli_abort("{.arg class_type} must select at least one class type.")
  }

  # Step 2: construct class-regression matrices in the exact unit order used
  # by the generated subject and marker random-effect arrays. The
  # longitudinal class-design scaffold crosses every subject with every
  # marker. This is important: an subject covariate such as treatment is
  # not marker-constant merely because the one-row model-matrix prototype came
  # from a single subject.
  marker_levels <- levels(
    factor(marker_prototype[[marker_variable]])
  ) # marker allocation units in their established factor order
  n_markers <- length(marker_levels) # number of marker allocation units
  class_long_data <- data_event[
    rep(seq_len(nrow(data_event)), each = n_markers),
    ,
    drop = FALSE
  ] # subject covariates repeated over every marker
  class_long_data[[marker_variable]] <- factor(
    rep(marker_levels, times = nrow(data_event)),
    levels = marker_levels
  ) # marker label for every subject-by-marker scaffold row
  prototype_only_columns <- setdiff(
    names(marker_prototype),
    names(class_long_data)
  ) # design variables defined only by the longitudinal prototype
  for (variable in prototype_only_columns) {
    marker_values <- marker_prototype[[variable]]
    if (length(marker_values) == n_markers) {
      class_long_data[[variable]] <- rep(
        marker_values,
        times = nrow(data_event)
      )
    }
  }
  class_design <- .mixture_class_design(
    formula_class = specification$formulaClass,
    selected_types = selected_types,
    data_long = class_long_data,
    data_event = data_event,
    id_variable = id_variable,
    marker_variable = marker_variable,
    n_classes = n_classes
  ) # subject- and marker-domain class design records
  ordering <- .resolve_class_ordering(
    class_ordering = specification$class_ordering,
    class_specific_formulae = isTRUE(class_design$class_specific)
  ) # contextual or explicit label-identification convention

  # Step 3: reuse the fitting data builder for coordinate validation, packed
  # starts, default coordinate choices and the shared-allocation description.
  fit_layout <- .build_mixture_standata(
    stan_data = list(
      n_id = nrow(data_event),
      D = n_markers,
      R_id = as.integer(dimensions[["subject"]]),
      R_mk = as.integer(dimensions[["marker"]]),
      Q_idm = as.integer(dimensions[["covariance_basis"]]),
      indep_idmarker_cov = as.integer(covariance_is_diagonal),
      zid_cols = labels$subject,
      zmk_cols = labels$marker,
      zidm_cols = labels$covariance
    ),
    mixture = list(
      n_classes = n_classes,
      class_type = selected_types,
      class_dimensions = specification$class_dimensions,
      class_prior = list(
        baseline_prob = specification$class_prior$baseline_prob,
        slope = specification$class_prior$slope,
        family = specification$class_prior$family
      ),
      class_design = class_design,
      class_ordering = ordering,
      include_survival = TRUE
    )
  ) # fit-identical coordinate layout and active-domain metadata

  # Step 4: resolve baseline probabilities. Probability ordering is a
  # constraint on supplied truth rather than a post-hoc relabelling step,
  # because relabelling class-specific regression formulae would change their
  # scientific meaning.
  class_parameters <- specification$class_parameters %||% list()
  if (!is.list(class_parameters)) {
    cli::cli_abort("{.arg class_parameters} must be a named list.")
  }
  if (
    length(class_parameters) > 0L &&
      (
        is.null(names(class_parameters)) ||
          any(!nzchar(names(class_parameters)))
      )
  ) {
    cli::cli_abort("{.arg class_parameters} must be a named list.")
  }
  allowed_parameter_names <- c(
    "probability",
    "probabilities",
    "location",
    "locations",
    "scale",
    "scales"
  )
  unknown_parameter_names <- setdiff(
    names(class_parameters) %||% character(0),
    allowed_parameter_names
  )
  if (length(unknown_parameter_names) > 0L) {
    cli::cli_abort(
      "Unknown {.arg class_parameters} entr{?y/ies}: {paste(unknown_parameter_names, collapse = ', ')}."
    )
  }
  alias_pairs <- list(
    c("probability", "probabilities"),
    c("location", "locations"),
    c("scale", "scales")
  ) # singular and plural convenience aliases which must not be duplicated
  duplicated_aliases <- vapply(alias_pairs, function(pair) {
    all(pair %in% names(class_parameters))
  }, logical(1))
  if (any(duplicated_aliases)) {
    offending_pair <- alias_pairs[[which(duplicated_aliases)[1L]]]
    cli::cli_abort(
      "Supply only one of {.arg {offending_pair[[1L]]}} and {.arg {offending_pair[[2L]]}} in {.arg class_parameters}."
    )
  }
  probability_supplied <- !is.null(
    class_parameters$probability %||% class_parameters$probabilities
  ) # whether strict probability order must be checked against user input
  probability_input <- class_parameters$probability %||%
    class_parameters$probabilities
  probability <- if (is.null(probability_input)) {
    drawn_probability <- stats::rgamma(n_classes, shape = 2, rate = 2)
    if (identical(ordering, "probability")) sort(drawn_probability) else drawn_probability
  } else {
    as.numeric(probability_input)
  } # unnormalised zero-covariate component probabilities, drawn once when omitted
  if (
    length(probability) != n_classes ||
      any(!is.finite(probability)) ||
      any(probability <= 0)
  ) {
    cli::cli_abort(
      "{.arg class_parameters$probability} must contain {n_classes} positive finite values."
    )
  }
  probability <- probability / sum(probability) # probability simplex used for sampling
  if (identical(ordering, "probability") && any(diff(probability) <= 0)) {
    detail <- if (probability_supplied) {
      "Supply strictly increasing probabilities in class-label order."
    } else {
      "The internally generated probabilities must be strictly increasing."
    }
    cli::cli_abort(c(
      x = "{.arg class_ordering = \"probability\"} requires strictly increasing baseline probabilities.",
      i = detail
    ))
  }

  # Step 5: align class-regression coefficients and calculate unit-specific
  # probabilities with the same stable multinomial-logit evaluator used for
  # posterior dynamic prediction.
  subject_coefficient_names <- class_design$subject$columns %||% character(0) # formulaClass coefficients governing subject allocations
  marker_coefficient_names <- class_design$marker$columns %||% character(0) # formulaClass coefficients governing marker allocations
  active_coefficient_domains <- c(
    subject = length(subject_coefficient_names) > 0L,
    marker = length(marker_coefficient_names) > 0L
  ) # domains contributing class-regression coefficients
  complete_coefficient_names <- if (sum(active_coefficient_domains) > 1L) {
    c(
      paste0("subject::", subject_coefficient_names),
      paste0("marker::", marker_coefficient_names)
    )
  } else {
    c(subject_coefficient_names, marker_coefficient_names)
  } # unambiguous coefficient order shared by simulation and fitting
  complete_class_coefficient <- .sim_resolve_prior_or_fixed(
    specification$class_slope_declaration,
    complete_coefficient_names,
    "class$slope"
  ) # one population draw or one fixed vector for the complete formulaClass regression
  number_subject_coefficients <- length(subject_coefficient_names) # split point between subject and marker domains
  coefficient_by_domain <- list(
    subject = stats::setNames(
      complete_class_coefficient[seq_len(number_subject_coefficients)],
      subject_coefficient_names
    ),
    marker = stats::setNames(
      complete_class_coefficient[number_subject_coefficients + seq_along(marker_coefficient_names)],
      marker_coefficient_names
    )
  ) # domain-specific views used by the multinomial-logit probability evaluator
  probability_by_domain <- lapply(c("subject", "marker"), function(domain) {
    probability_array <- .mixture_class_probability(
      baseline_probability = matrix(probability, nrow = 1L),
      class_coefficient = matrix(
        coefficient_by_domain[[domain]],
        nrow = 1L
      ),
      class_design = class_design[[domain]]$matrix,
      class_term_start = class_design[[domain]]$class_term_start,
      class_term_count = class_design[[domain]]$class_term_count
    ) # one-draw unit-by-class probabilities
    matrix(
      probability_array[1L, , ],
      nrow = nrow(class_design[[domain]]$matrix),
      ncol = n_classes,
      dimnames = list(
        class_design[[domain]]$units,
        paste0("class_", seq_len(n_classes))
      )
    )
  })
  names(probability_by_domain) <- c("subject", "marker")

  # Step 6: draw exactly one class for each active natural allocation domain.
  # Compatible blocks refer to the same vector, preventing accidental
  # Cartesian-product allocations.
  allocation_domains <- fit_layout$mixture$allocation_domains
  allocation <- list(
    subject = if (length(allocation_domains$subject) > 0L) {
      stats::setNames(
        vapply(seq_len(nrow(probability_by_domain$subject)), function(unit) {
          sample.int(
            n_classes,
            size = 1L,
            prob = probability_by_domain$subject[unit, ]
          )
        }, integer(1)),
        rownames(probability_by_domain$subject)
      )
    } else {
      integer(0)
    },
    marker = if (length(allocation_domains$marker) > 0L) {
      stats::setNames(
        vapply(seq_len(nrow(probability_by_domain$marker)), function(unit) {
          sample.int(
            n_classes,
            size = 1L,
            prob = probability_by_domain$marker[unit, ]
          )
        }, integer(1)),
        rownames(probability_by_domain$marker)
      )
    } else {
      integer(0)
    }
  ) # realised latent classes, never appended to the observed datasets

  # Step 7: resolve complete packed location and scale matrices. Omitted
  # population matrices are drawn once, then fixed for all allocation units.
  total_dimension <- fit_layout$K_mix # packed selected-coordinate count
  default_centre <- seq(-1.25, 1.25, length.out = n_classes)
  default_location <- matrix(
    rep(default_centre, total_dimension),
    nrow = n_classes,
    ncol = total_dimension
  ) + matrix(
    stats::rnorm(n_classes * total_dimension, mean = 0, sd = 0.15),
    nrow = n_classes,
    ncol = total_dimension
  ) # population component locations drawn once around separated class centres
  ordered_coordinate <- fit_layout$mix_ordered_location_coordinate
  if (identical(ordering, "intercept") && ordered_coordinate > 0L) {
    default_location[, ordered_coordinate] <- sort(default_location[, ordered_coordinate])
  }
  location <- .sim_mixture_component_matrix(
    class_parameters$location %||% class_parameters$locations,
    default = default_location,
    layout = fit_layout,
    argument_name = "class_parameters$location"
  ) # complete component-location matrix
  scale <- .sim_mixture_component_matrix(
    class_parameters$scale %||% class_parameters$scales,
    default = matrix(
      stats::runif(n_classes * total_dimension, min = 0.45, max = 0.85),
      nrow = n_classes,
      ncol = total_dimension
    ),
    layout = fit_layout,
    argument_name = "class_parameters$scale",
    require_positive = TRUE
  ) # complete positive component-scale matrix
  if (
    identical(ordering, "intercept") &&
      ordered_coordinate > 0L &&
      any(diff(location[, ordered_coordinate]) <= 0)
  ) {
    cli::cli_abort(c(
      x = "The ordered class-location coordinate is not strictly increasing.",
      i = "Arrange {.arg class_parameters$location} in class-label order."
    ))
  }

  list(
    n_classes = n_classes,
    class_type = fit_layout$mixture$class_type,
    dimensions = fit_layout$mixture$dimensions,
    starts = fit_layout$mixture$starts,
    total_dimension = total_dimension,
    ordering = ordering,
    ordered_location_coordinate = ordered_coordinate,
    distribution = fit_layout$mixture$distribution,
    component_family = specification$class_prior$family,
    probability = stats::setNames(
      probability,
      paste0("class_", seq_len(n_classes))
    ),
    coefficient = coefficient_by_domain,
    probability_by_unit = probability_by_domain,
    allocation = allocation,
    location = location,
    scale = scale,
    class_design = class_design,
    allocation_domains = allocation_domains,
    stan_data_contract = list(
      use_mixture = fit_layout$use_mixture,
      n_classes = fit_layout$n_classes,
      K_mix = fit_layout$K_mix,
      mix_subject = fit_layout$mix_subject,
      mix_marker = fit_layout$mix_marker,
      mix_covariance = fit_layout$mix_covariance,
      mix_dim_subject = fit_layout$mix_dim_subject,
      mix_idx_subject = fit_layout$mix_idx_subject,
      mix_start_subject = fit_layout$mix_start_subject,
      mix_dim_marker = fit_layout$mix_dim_marker,
      mix_idx_marker = fit_layout$mix_idx_marker,
      mix_start_marker = fit_layout$mix_start_marker,
      mix_dim_covariance = fit_layout$mix_dim_covariance,
      mix_idx_covariance = fit_layout$mix_idx_covariance,
      mix_start_covariance = fit_layout$mix_start_covariance,
      mix_ordering = fit_layout$mix_ordering,
      mix_ordered_location_coordinate =
        fit_layout$mix_ordered_location_coordinate,
      P_class_subject = fit_layout$P_class_subject,
      X_class_subject = fit_layout$X_class_subject,
      class_term_start_subject =
        fit_layout$class_term_start_subject,
      class_term_count_subject =
        fit_layout$class_term_count_subject,
      P_class_marker = fit_layout$P_class_marker,
      X_class_marker = fit_layout$X_class_marker,
      class_term_start_marker =
        fit_layout$class_term_start_marker,
      class_term_count_marker =
        fit_layout$class_term_count_marker,
      mix_component_family = fit_layout$mix_component_family,
      mix_component_df = fit_layout$mix_component_df
    ) # structural Stan fields which a fit to the simulated data must reproduce exactly
  )
}

#' Apply a simulated component distribution to selected coordinates
#'
#' @param latent_matrix Allocation-unit by coordinate matrix initially drawn
#'   from the ordinary standardised distribution.
#' @param level Canonical class-specific block name.
#' @param allocation Integer class label per row.
#' @param mixture Prepared generative mixture, or `NULL`.
#'
#' @return The input matrix with selected coordinates replaced by conditional
#'   component draws.
#' @keywords internal
#' @noRd
.sim_apply_mixture_to_latent <- function(
  latent_matrix,
  level,
  allocation,
  mixture
) {
  if (is.null(mixture) || !(level %in% mixture$class_type)) {
    return(latent_matrix)
  }
  selected_coordinates <- mixture$dimensions[[level]] # source coordinates to replace
  packed_start <- mixture$starts[[level]] # first corresponding mixture column
  if (length(selected_coordinates) == 0L) {
    return(latent_matrix)
  }
  if (nrow(latent_matrix) != length(allocation)) {
    cli::cli_abort(
      "Internal mixture allocation does not align with the {level} latent matrix."
    )
  }

  # Each coordinate is sampled conditionally and independently, matching the
  # product component density used in Stan. The allocation itself remains
  # common across all selected coordinates and compatible random-effect blocks.
  for (local_coordinate in seq_along(selected_coordinates)) {
    source_coordinate <- selected_coordinates[[local_coordinate]]
    packed_coordinate <- packed_start + local_coordinate - 1L
    standard_draw <- .sim_draw_standard_component(
      nrow(latent_matrix),
      family = mixture$component_family
    ) # centred unit-scale draw from the class family declared through jm_truth() and jm_priors()
    latent_matrix[, source_coordinate] <-
      mixture$location[allocation, packed_coordinate] +
      mixture$scale[allocation, packed_coordinate] * standard_draw
  }
  latent_matrix
}
