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
#'
#' @details
#' Let \(C_j\in\{1,\ldots,G\}\) denote the latent class for allocation unit
#' \(j\). Conditional on class \(g\), every selected coordinate \(r\) is drawn
#' independently as
#'
#' \deqn{
#' z_{jr}\mid C_j=g \sim
#' D\{\mu_{gr},\sigma_{gr}\},
#' }
#'
#' where \(D\) is Student-\(t_6\), Laplace or Normal according to `shrinkage`.
#' Unselected subject, marker and covariance-regression coordinates remain
#' standard Normal.
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
#'   to one. The default is uniform, except under probability ordering where a
#'   strictly increasing sequence is used.
#' - `coefficient`: class-regression coefficients. Supply a named list with
#'   `subject` and/or `marker` vectors. Named vectors are matched to the
#'   columns recorded in `truth$mixture$class_design`; unnamed vectors must
#'   have the exact required length.
#' - `location`: an `n_classes` by \(K\) matrix in the packed coordinate
#'   order, or a named list containing matrices for selected class
#'   types. The default places ordered, equally spaced component centres
#'   between -1.25 and 1.25 on every selected coordinate.
#' - `scale`: a positive scalar, vector, matrix, or named level list matching
#'   `location`. The default component scale is 0.65.
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
#' @param class_parameters Generative class probabilities, regression
#'   coefficients, component locations and component scales. See Details.
#'
#' @return A list with the same top-level components as [simulate_joinme()].
#'   The `truth` and `true_params` entries additionally contain a `mixture`
#'   record suitable for recovery studies and class-allocation assessment.
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
#'   class_parameters = list(
#'     probability = c(0.25, 0.45, 0.30),
#'     coefficient = list(subject = c(
#'       "class_1:x1" = -0.6,
#'       "class_2:x1" = 0.4
#'     )),
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
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  ...,
  seed = .Random.seed[[1]],
  covariate_formulas = list(
    x1 ~ rnorm(n_id),
    x2 ~ rnorm(n_id)
  ),
  marker_weights = NULL,
  shared_marker_weights = TRUE,
  fixed_marker_weights = FALSE,
  shrinkage = 0L,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6),
  beta_long = NULL,
  beta_event = NULL,
  dist_coefs = list(),
  re_params = list(
    id = list(sd = NULL, corr = NULL),
    marker = list(sd = NULL, corr = NULL),
    id_marker_cov = list(
      latent = list(sd = NULL, corr = NULL),
      alpha = NULL,
      beta = NULL,
      lambda = NULL,
      diag_link = "softplus"
    ),
    dist = list()
  ),
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
  h0 = NULL,
  baseline_hazard = list(type = "weibull", shape = 1.4, scale = 6.0),
  formulaBasehaz = NULL,
  beta_basehaz = NULL,
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
  selected_class_types <- .canonical_class_types(
    class_type
  ) # checked public class types before ordinary simulation work begins
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

  # Step 1: retain the fitting entry point's mixture syntax in one small,
  # unevaluated specification. The ordinary simulator will validate the
  # formula-derived dimensions after it has built the same random-effect
  # design matrices used for data generation.
  mixture_specification <- list(
    n_classes = n_classes, # requested number of latent components
    formulaClass = formulaClass, # covariate model for component membership
    class_type = selected_class_types, # random-effect blocks receiving mixture distributions
    class_dimensions = class_dimensions, # selected coordinates by block
    class_ordering = class_ordering, # label-identification convention
    class_parameters = class_parameters # fixed generative mixture parameters
  )

  # Step 2: reconstruct the ordinary simulation call from the arguments that
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

  # A longitudinal-only request still needs subject covariates while the
  # ordinary engine constructs model matrices. A zero-hazard internal
  # survival formula provides that scaffold without censoring longitudinal
  # observations or allowing the event process to affect any generated value.
  if (longitudinal_only) {
    simulation_arguments$formulaEvent <-
      survival::Surv(time, event) ~ 1
    simulation_arguments$formulaAssoc <- NULL
    simulation_arguments$assoc <- "cv_total"
    simulation_arguments$assoc_coefs <- c(cv_total = 0)
    simulation_arguments$h0 <- function(t) {
      rep(0, length(t))
    }
    simulation_arguments$censor_longitudinal_after_event <- FALSE
    simulation_arguments$root_control <- list(
      t_init = time_cens,
      t_max = time_cens,
      expand = 1.7,
      max_expand = 1L
    )
  }

  # Step 3: execute the common generator. The private specification changes
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
    simulation$true_params$formulaEvent <- NULL
  }
  simulation
}

#' Resolve one domain's class-regression coefficients for simulation
#'
#' @param supplied User-supplied numeric vector, or `NULL`.
#' @param design_record Allocation-domain design produced by
#'   `.mixture_class_design()`.
#' @param domain Plain-language allocation-domain name.
#'
#' @return A named numeric vector in the fitted concatenated-column order.
#' @keywords internal
#' @noRd
.sim_mixture_class_coefficient <- function(supplied, design_record, domain) {
  required_names <- design_record$columns %||% character(0) # fitted coefficient labels
  required_count <- length(required_names) # number of coefficients for this domain
  if (required_count == 0L) {
    if (!is.null(supplied) && length(supplied) > 0L) {
      cli::cli_abort(
        "No class-regression coefficient is required for the {domain} domain."
      )
    }
    return(stats::setNames(numeric(0), character(0)))
  }
  if (is.null(supplied)) {
    return(stats::setNames(rep(0, required_count), required_names))
  }

  supplied_names <- names(supplied) # optional design-column labels
  supplied <- as.numeric(supplied) # finite coefficient values to be aligned
  if (any(!is.finite(supplied))) {
    cli::cli_abort(
      "{.arg class_parameters$coefficient} for {domain} must be finite."
    )
  }
  if (!is.null(supplied_names) && all(nzchar(supplied_names))) {
    unknown_names <- setdiff(supplied_names, required_names)
    if (length(unknown_names) > 0L) {
      cli::cli_abort(
        "Unknown {domain} class coefficient{?s}: {paste(unknown_names, collapse = ', ')}."
      )
    }
    coefficient <- stats::setNames(rep(0, required_count), required_names)
    coefficient[supplied_names] <- supplied
    return(coefficient)
  }
  if (length(supplied) != required_count) {
    cli::cli_abort(c(
      x = "{.arg class_parameters$coefficient} has the wrong length for {domain}.",
      i = "Expected {required_count} value{?s}; received {length(supplied)}."
    ))
  }
  stats::setNames(supplied, required_names)
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
#' @param shrinkage Integer component-family code.
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
  shrinkage,
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
      shrinkage = as.integer(shrinkage),
      zid_cols = labels$subject,
      zmk_cols = labels$marker,
      zidm_cols = labels$covariance
    ),
    mixture = list(
      n_classes = n_classes,
      class_type = selected_types,
      class_dimensions = specification$class_dimensions,
      class_probability_concentration = 1,
      class_regression_scale = 1,
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
    "coefficient",
    "coefficients",
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
    c("coefficient", "coefficients"),
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
  probability <- as.numeric(
    class_parameters$probability %||%
      class_parameters$probabilities %||%
      if (identical(ordering, "probability")) seq_len(n_classes) else rep(1, n_classes)
  ) # unnormalised zero-covariate component probabilities
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
  coefficient_input <-
    class_parameters$coefficient %||% class_parameters$coefficients
  if (
    is.list(coefficient_input) &&
      (
        is.null(names(coefficient_input)) ||
          any(!nzchar(names(coefficient_input)))
      )
  ) {
    cli::cli_abort(
      "{.arg class_parameters$coefficient} must use subject and/or marker names."
    )
  }
  if (is.list(coefficient_input)) {
    unknown_coefficient_domains <- setdiff(
      names(coefficient_input),
      c("subject", "marker")
    )
    if (length(unknown_coefficient_domains) > 0L) {
      cli::cli_abort(
        "Unknown class-coefficient domain{?s}: {paste(unknown_coefficient_domains, collapse = ', ')}."
      )
    }
  }
  if (
    !is.null(coefficient_input) &&
      !is.list(coefficient_input) &&
      sum(c(
        ncol(class_design$subject$matrix) > 0L,
        ncol(class_design$marker$matrix) > 0L
      )) > 1L
  ) {
    cli::cli_abort(c(
      x = "{.arg class_parameters$coefficient} is ambiguous across two allocation domains.",
      i = "Supply list(subject = ..., marker = ...)."
    ))
  }
  coefficient_by_domain <- list(
    subject = .sim_mixture_class_coefficient(
      if (is.list(coefficient_input)) coefficient_input$subject else coefficient_input,
      class_design$subject,
      "subject"
    ),
    marker = .sim_mixture_class_coefficient(
      if (is.list(coefficient_input)) coefficient_input$marker else coefficient_input,
      class_design$marker,
      "marker"
    )
  ) # concatenated regression coefficients in design-column order
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

  # Step 7: resolve complete packed location and scale matrices. Defaults are
  # intentionally separated enough for simulation examples while remaining
  # moderate on the standardised random-effect scale.
  total_dimension <- fit_layout$K_mix # packed selected-coordinate count
  default_centre <- seq(-1.25, 1.25, length.out = n_classes)
  default_location <- matrix(
    rep(default_centre, total_dimension),
    nrow = n_classes,
    ncol = total_dimension
  ) # default ordered component centres on every selected coordinate
  location <- .sim_mixture_component_matrix(
    class_parameters$location %||% class_parameters$locations,
    default = default_location,
    layout = fit_layout,
    argument_name = "class_parameters$location"
  ) # complete component-location matrix
  scale <- .sim_mixture_component_matrix(
    class_parameters$scale %||% class_parameters$scales,
    default = matrix(0.65, nrow = n_classes, ncol = total_dimension),
    layout = fit_layout,
    argument_name = "class_parameters$scale",
    require_positive = TRUE
  ) # complete positive component-scale matrix
  ordered_coordinate <- fit_layout$mix_ordered_location_coordinate
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
      shrinkage = as.integer(shrinkage)
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
#' @param shrinkage Integer component-family code.
#'
#' @return The input matrix with selected coordinates replaced by conditional
#'   component draws.
#' @keywords internal
#' @noRd
.sim_apply_mixture_to_latent <- function(
  latent_matrix,
  level,
  allocation,
  mixture,
  shrinkage
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
    standard_draw <- .sim_draw_standard_shrinkage(
      nrow(latent_matrix),
      shrinkage = shrinkage
    ) # centred unit-scale draw from the selected component family
    latent_matrix[, source_coordinate] <-
      mixture$location[allocation, packed_coordinate] +
      mixture$scale[allocation, packed_coordinate] * standard_draw
  }
  latent_matrix
}
