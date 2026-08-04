#' Latent-progress mixtures for nested multivariate joint models
#'
#' @description
#' `joinme_mix()` fits the JoiNMe model with a finite mixture on one or more
#' standardised random-effect blocks.  It has a separate entry point and
#' fitted-object class, whilst deliberately retaining the formula parser,
#' likelihood, association transforms, prediction programme, diagnostics, and
#' posterior utilities used by [joinme()].
#'
#' The allocation variables are analytically marginalised in Stan.  This is
#' important because Hamiltonian Monte Carlo cannot sample discrete class
#' labels.  Posterior class probabilities are reconstructed for every
#' subject or marker in generated quantities.
#'
#' @details
#' The longitudinal linear predictor remains
#'
#' \deqn{
#' \eta_{id}(t)=x_{id}(t)^\top\beta+
#' z_i(t)^\top u_i+z_d(t)^\top v_d+
#' z_{id}(t)^\top w_{id}.
#' }
#'
#' Latent-class modelling changes the distribution of selected *standardised* random
#' effects, not the structural decomposition above.  For a selected vector
#' \eqn{a_j}, the ordinary density \eqn{D(a_j;0,1)} is replaced by
#'
#' \deqn{
#' p(a_j)=\sum_{g=1}^{G}\pi_g
#' D\{a_j;\mu_g,\operatorname{diag}(s_g)\}.
#' }
#'
#' `shrinkage = 0`, `1`, and `2` select Student-\eqn{t_6}, Laplace, and Normal
#' latent-class component densities, respectively. These choices do not set
#' the ordinary priors for `beta`, `alpha`, `iota`, marker effects or marker
#' weights; declare those independently with [jm_prior()]. The component locations and scales are
#' estimated. By default, Stan orders only the first selected random-intercept
#' location. Selected slopes and other coordinates remain unrestricted.
#' Alternatively, baseline class probabilities may be ordered, or ordering may
#' be disabled with `class_ordering`.
#'
#' The supported class types are:
#'
#' - `"subject"`: the standardised subject effects underlying \eqn{u_i};
#' - `"marker"`: the standardised marker effects underlying \eqn{v_d};
#' - `"corr"`: the off-diagonal correlation coordinates in the
#'   marker-by-subject covariance regression;
#' - `"vcov"`: the lower-triangular variance--covariance coordinates in the
#'   marker-by-subject covariance regression.
#'
#' Compatible types share a single class allocation. In particular,
#' `c("subject", "vcov")` gives one \eqn{G}-class probability vector per subject.
#' The construction therefore has \eqn{G}
#' classes, rather than the Cartesian product of separate level-specific
#' classes. If both subject-indexed and marker-indexed types are selected,
#' they use the same \eqn{G} component labels and mixing proportions but retain
#' the natural allocation unit of their block.
#'
#' By default the first two coordinates of each random-effect block define the
#' two-dimensional progress plane.  A one-dimensional block is retained as a
#' well-defined special case.  `class_dimensions` can select other
#' coordinates. `corr` and `vcov` are alternatives because both act on the
#' same covariance-regression block.
#'
#' If `formulaEvent` and `dataEvent` are both `NULL`, a likelihood-neutral event
#' scaffold is created solely to reuse the common design code.  Event times and
#' event indicators are then set to zero in Stan data, so inference is based
#' only on the longitudinal process. Association terms are consequently
#' unavailable in this mode.
#'
#' @param formulaLong,dataLong,formulaVCov,formulaDist,control,draws,families,
#'   transforms,priors,fixed_marker_weights,shared_marker_weights,basehaz,fit,
#'   seed Arguments with the same meaning as in [joinme()].
#' @param formulaEvent Optional survival formula.  Supply this together with
#'   `dataEvent`, or leave both `NULL` for a longitudinal-only mixture.
#' @param dataEvent Optional event-process data.
#' @param n_classes Number of latent classes \eqn{G}; must be at least two.
#' @param class_type Character vector naming the random-effect class
#'   types. The only accepted values are `"subject"`, `"marker"`, `"corr"`,
#'   and `"vcov"`.
#' @param class_dimensions Optional named list of integer coordinate indices,
#'   one entry per selected class type. For a single type, an integer vector is
#'   also accepted directly.
#' @param formulaClass A one-sided formula for class-membership covariates, or a
#'   named list with `subject` and `marker` formulae. A list containing
#'   exactly `n_classes` formulae gives a separate predictor for every class.
#'   Covariates must be constant within the corresponding allocation unit. The
#'   intercept is represented by the baseline class probabilities.
#' @param class_ordering Identification rule for the common class labels.
#'   `"intercept"` orders only the first selected random-intercept location;
#'   `"probability"` orders the baseline class probabilities; `"none"` imposes
#'   no order. The default is `"intercept"` for a shared class formula and
#'   `"none"` for a list of `n_classes` class-specific formulae.
#' @param ... Additional arguments passed to [joinme_standata()], including
#'   `assoc`, `marker_weights`, `id_var`, `marker_var`, `time_var`, and
#'   `shrinkage`.
#'
#' @return If `fit = TRUE`, a `JoiNMeMixFit` object inheriting from
#'   `JoiNMeFit`.  If `fit = FALSE`, a prepared `JoiNMeMixStanData` object.
#'
#' @examples
#' \dontrun{
#' mixed_fit <- joinme_mix(
#'   y ~ time + (1 + time | id) +
#'     (1 + time + (1 + time | id) | marker),
#'   dataLong,
#'   survival::Surv(time, event) ~ treatment,
#'   dataEvent,
#'   n_classes = 3,
#'   class_type = c("subject", "vcov"),
#'   formulaClass = ~ treatment + age,
#'   priors = jm_prior(class_probability = rep(2, 3)),
#'   assoc = c("cv_total")
#' )
#'
#' longitudinal_only_fit <- joinme_mix(
#'   y ~ time + (1 + time | id) +
#'     (1 + time + (1 + time | id) | marker),
#'   dataLong,
#'   n_classes = 2,
#'   class_type = "subject"
#' )
#' }
#' @export
joinme_mix <- function(
  formulaLong,
  dataLong,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaVCov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = joinme_priors(),
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  basehaz = joinme_basehaz(),
  n_classes = 2L,
  formulaClass = ~1,
  class_type = "subject",
  class_dimensions = NULL,
  class_ordering = NULL,
  fit = TRUE,
  seed = NULL,
  ...
) {
  prior_specification <- unclass(
    .joinme_priors_(priors, validate = TRUE)
  ) # validated common and latent-class prior declarations

  # Step 1: distinguish a genuine survival analysis from a longitudinal-only
  # analysis.  Requiring the formula and data together prevents a partly
  # specified event process from being mistaken for a censored analysis.
  has_survival_process <- !is.null(formulaEvent) || !is.null(dataEvent)
  if (xor(is.null(formulaEvent), is.null(dataEvent))) {
    cli::cli_abort(c(
      x = "{.arg formulaEvent} and {.arg dataEvent} must be supplied together.",
      i = "Leave both NULL for a longitudinal-only latent-progress model."
    ))
  }

  dot_arguments <- list(...)
  association_terms <- dot_arguments$assoc
  if (is.null(association_terms)) {
    association_terms <- if (has_survival_process) "cv_mean" else character(0)
  }

  # Step 2: a longitudinal-only fit still needs the common matrix-building
  # machinery to know the study time scale.  The scaffold copies one observed
  # row per subject so every longitudinal covariate used in a design
  # template remains available.  Its survival contribution is removed later by
  # setting both integration limits to zero.
  if (!has_survival_process) {
    if (length(association_terms) > 0L) {
      cli::cli_abort(c(
        x = "Association terms require a survival process.",
        i = "Remove {.arg assoc}, or supply both {.arg formulaEvent} and {.arg dataEvent}."
      ))
    }
    event_scaffold <- .longitudinal_only_event_scaffold(
      data_long = dataLong,
      id_variable = dot_arguments$id_var %||% "id",
      time_variable = dot_arguments$time_var %||% "time"
    )
    formulaEvent <- survival::Surv(
      .joinme_follow_up,
      .joinme_event
    ) ~ 1
    dataEvent <- event_scaffold
  }

  # Step 3: retain the user's mixture request as a plain specification.  Its
  # dimensions can only be checked after the ordinary formula parser has built
  # the four random-effect blocks, so final validation occurs in
  # `.build_mixture_standata()`.
  if (
    !is.numeric(n_classes) ||
      length(n_classes) != 1L ||
      !is.finite(n_classes) ||
      n_classes < 2 ||
      n_classes != floor(n_classes)
  ) {
    cli::cli_abort("{.arg n_classes} must be a single integer of at least two.")
  }
  n_classes <- as.integer(n_classes) # validated number of latent classes
  selected_types <- .canonical_class_types(
    class_type
  ) # canonical random-effect blocks requested for latent-class modelling
  class_design <- .mixture_class_design(
    formula_class = formulaClass,
    selected_types = selected_types,
    data_long = dataLong,
    data_event = dataEvent,
    id_variable = dot_arguments$id_var %||% "id",
    marker_variable = dot_arguments$marker_var %||% "marker",
    n_classes = n_classes
  ) # allocation-unit design matrices for class-probability regression
  resolved_class_ordering <- .resolve_class_ordering(
    class_ordering = class_ordering,
    class_specific_formulae = isTRUE(class_design$class_specific)
  ) # explicit or context-sensitive class-label identification rule

  mixture_specification <- list(
    n_classes = n_classes,
    class_type = selected_types,
    class_dimensions = class_dimensions,
    class_probability_concentration =
      prior_specification$class_probability,
    class_regression_prior = prior_specification$class_regression,
    class_design = class_design,
    class_ordering = resolved_class_ordering,
    include_survival = has_survival_process
  )

  dot_arguments$assoc <- association_terms
  dot_arguments$mixture <- mixture_specification

  # Step 4: invoke the established JoiNMe preparation path.  `fit = FALSE`
  # deliberately stops immediately before sampling; the returned holder is
  # then marked so its existing sampler constructs the mixture subclass.
  joinme_arguments <- c(
    list(
      formulaLong = formulaLong,
      dataLong = dataLong,
      formulaEvent = formulaEvent,
      dataEvent = dataEvent,
      formulaVCov = formulaVCov,
      formulaDist = formulaDist,
      control = control,
      draws = draws,
      families = families,
      transforms = transforms,
      priors = priors,
      fixed_marker_weights = fixed_marker_weights,
      shared_marker_weights = shared_marker_weights,
      basehaz = basehaz,
      fit = FALSE,
      seed = seed
    ),
    dot_arguments
  )

  prepared_model <- do.call(joinme, joinme_arguments)
  result_metadata <- prepared_model$stan_data$mixture
  prepared_model$parent_call <- match.call()
  prepared_model <- JoiNMeMixStanData$new(
    prepared_model,
    more_recipe = list(mixture = result_metadata)
  )
  prepared_model$result_class <- "mixture" # explicit fitted-result contract stored on the final mixture preparation holder

  if (!isTRUE(fit)) {
    return(prepared_model)
  }
  prepared_model$sample()
}

#' Resolve latent-class label ordering
#'
#' @param class_ordering User-supplied ordering rule, or `NULL`.
#' @param class_specific_formulae Whether every class has its own formula.
#'
#' @return One of `"intercept"`, `"probability"`, or `"none"`.
#' @keywords internal
#' @noRd
.resolve_class_ordering <- function(
  class_ordering,
  class_specific_formulae
) {
  if (is.null(class_ordering)) {
    return(if (isTRUE(class_specific_formulae)) "none" else "intercept")
  }
  if (!is.character(class_ordering) || length(class_ordering) != 1L) {
    cli::cli_abort(
      "{.arg class_ordering} must be intercept, probability, or none."
    )
  }
  match.arg(
    tolower(class_ordering),
    c("intercept", "probability", "none")
  )
}

#' Construct allocation-unit class-regression designs
#'
#' @param formula_class One formula shared by active domains, or a named list of
#'   domain-specific formulae. A list of `n_classes` formulae gives one formula
#'   per class.
#' @param n_classes Number of latent classes.
#' @param selected_types Canonical random-effect levels assigned latent classes.
#' @param data_long Longitudinal data supplied to `joinme_mix()`.
#' @param data_event Event data, including the longitudinal-only scaffold.
#' @param id_variable Subject identifier column.
#' @param marker_variable Marker identifier column.
#'
#' @return A named list containing formulae, matrices, unit labels and columns.
#' @keywords internal
#' @noRd
.mixture_class_design <- function(
  formula_class,
  selected_types,
  data_long,
  data_event,
  id_variable,
  marker_variable,
  n_classes = 2L
) {
  active_domains <- c(
    subject = any(selected_types %in% c("subject", "corr", "vcov")),
    marker = "marker" %in% selected_types
  ) # whether each natural allocation domain contributes to the mixture

  n_classes <- as.integer(n_classes) # number of class-specific formula positions
  is_group_formula_list <- function(value) {
    is.list(value) &&
      length(value) == n_classes &&
      all(vapply(value, inherits, logical(1), what = "formula"))
  }
  normalise_formula_specification <- function(value, domain) {
    if (inherits(value, "formula")) {
      return(list(
        formulae = c(
          rep(list(value), max(0L, n_classes - 1L)),
          list(~1)
        ),
        class_specific = FALSE
      ))
    }
    if (is_group_formula_list(value)) {
      return(list(
        formulae = unname(value),
        class_specific = TRUE
      ))
    }
    cli::cli_abort(c(
      x = "Invalid {.arg formulaClass} for the {domain} domain.",
      i = "Use one formula, or a list containing exactly {n_classes} one-sided formulae."
    ))
  }

  is_domain_formula_list <- is.list(formula_class) &&
    !is.null(names(formula_class)) &&
    all(names(formula_class) %in% names(active_domains))
  if (inherits(formula_class, "formula") ||
      (is_group_formula_list(formula_class) && !is_domain_formula_list)) {
    formulae <- list(
      subject = normalise_formula_specification(
        formula_class,
        "subject"
      ),
      marker = normalise_formula_specification(
        formula_class,
        "marker"
      )
    ) # one specification is shared when no domain names are supplied
  } else if (is_domain_formula_list) {
    unknown_domains <- setdiff(
      names(formula_class),
      names(active_domains)
    ) # formula-list entries which do not identify an allocation domain
    if (length(unknown_domains) > 0L) {
      cli::cli_abort(
        "Unknown {.arg formulaClass} domain(s): {paste(unknown_domains, collapse = ', ')}."
      )
    }
    formulae <- list(
      subject = normalise_formula_specification(
        formula_class$subject %||% ~1,
        "subject"
      ),
      marker = normalise_formula_specification(
        formula_class$marker %||% ~1,
        "marker"
      )
    ) # omitted domains retain constant baseline probabilities
  } else {
    cli::cli_abort(
      "{.arg formulaClass} must be a formula, a list of {n_classes} formulae, or a named domain list."
    )
  }

  build_domain <- function(
    formula_specification,
    primary_data,
    fallback_data,
    unit_variable,
    unit_levels,
    domain
  ) {
    class_formulae <- formula_specification$formulae
    valid_formulae <- vapply(
      class_formulae,
      function(formula) inherits(formula, "formula") && length(formula) == 2L,
      logical(1)
    ) # whether every class has a one-sided formula
    if (!all(valid_formulae)) {
      cli::cli_abort(
        "Every {.arg formulaClass} entry for the {domain} domain must be one-sided."
      )
    }
    variables <- unique(unlist(
      lapply(class_formulae, all.vars),
      use.names = FALSE
    )) # covariates requested by at least one class formula
    unit_frame <- data.frame(
      .class_unit = unit_levels,
      stringsAsFactors = FALSE
    ) # one row per allocation unit in Stan's established ordering

    for (variable in variables) {
      source_data <- if (variable %in% names(primary_data)) {
        primary_data
      } else if (
        !is.null(fallback_data) &&
          variable %in% names(fallback_data)
      ) {
        fallback_data
      } else {
        cli::cli_abort(
          "Class covariate {.val {variable}} is unavailable for the {domain} domain."
        )
      } # data source containing the requested baseline covariate

      values <- lapply(unit_levels, function(unit) {
        unit_rows <- source_data[[unit_variable]] == unit
        observed <- unique(source_data[[variable]][unit_rows])
        observed <- observed[!is.na(observed)]
        if (length(observed) != 1L) {
          cli::cli_abort(c(
            x = "Class covariate {.val {variable}} is not constant within {domain} unit {.val {unit}}.",
            i = "Class-regression covariates must be measured once, or remain constant within their allocation unit."
          ))
        }
        if (is.factor(observed)) {
          as.character(observed[[1L]])
        } else {
          observed[[1L]]
        }
      }) # checked unit-level covariate values
      unit_frame[[variable]] <- unlist(
        values,
        recursive = FALSE,
        use.names = FALSE
      )
      if (is.factor(source_data[[variable]])) {
        unit_frame[[variable]] <- factor(
          unit_frame[[variable]],
          levels = levels(source_data[[variable]])
        )
      }
    }

    formula_records <- lapply(seq_len(n_classes), function(group) {
      formula <- class_formulae[[group]]
      model_frame <- stats::model.frame(
        formula,
        data = unit_frame,
        na.action = stats::na.pass
      ) # allocation-unit data after applying this class formula
      formula_terms <- stats::terms(
        model_frame
      ) # fitted transformation record for this class
      factor_levels <- stats::.getXlevels(
        formula_terms,
        model_frame
      ) # categorical levels retained for prediction
      design <- stats::model.matrix(
        formula_terms,
        data = model_frame
      ) # conventional class-specific design
      contrast_rules <- attr(
        design,
        "contrasts"
      ) # fitted contrast rules for this class
      design <- design[
        ,
        setdiff(colnames(design), "(Intercept)"),
        drop = FALSE
      ] # the baseline probability supplies the class intercept
      list(
        class = group,
        formula = formula,
        terms = formula_terms,
        factor_levels = factor_levels,
        contrasts = contrast_rules,
        matrix = unname(design),
        columns = colnames(design)
      )
    }) # fitted design record for every latent class
    class_term_count <- vapply(
      formula_records,
      function(record) ncol(record$matrix),
      integer(1)
    ) # number of regression coefficients used by each class
    class_term_start <- cumsum(c(
      1L,
      head(class_term_count, -1L)
    )) # first concatenated coefficient position for each class
    concatenated_design <- do.call(
      cbind,
      lapply(formula_records, `[[`, "matrix")
    ) # compact unit design containing only coefficients actually requested
    if (is.null(concatenated_design)) {
      concatenated_design <- matrix(0, nrow = length(unit_levels), ncol = 0L)
    }
    class_term_start[class_term_count == 0L] <- pmin(
      class_term_start[class_term_count == 0L],
      max(1L, ncol(concatenated_design))
    ) # valid inert start positions for classes without coefficients
    coefficient_class <- rep(
      seq_len(n_classes),
      times = class_term_count
    ) # class receiving each concatenated coefficient
    coefficient_covariate <- unlist(
      lapply(formula_records, `[[`, "columns"),
      use.names = FALSE
    ) # fitted design-column label for each concatenated coefficient

    list(
      formulae = class_formulae,
      formula_records = formula_records,
      class_specific = isTRUE(formula_specification$class_specific),
      matrix = unname(concatenated_design),
      columns = paste0(
        "class_",
        coefficient_class,
        ":",
        coefficient_covariate
      ),
      coefficient_class = coefficient_class,
      coefficient_covariate = coefficient_covariate,
      class_term_start = as.integer(class_term_start),
      class_term_count = as.integer(class_term_count),
      units = as.character(unit_levels)
    )
  }

  subject_levels <- sort(
    unique(data_event[[id_variable]])
  ) # subject order used by joinme_standata()
  marker_levels <- levels(
    factor(data_long[[marker_variable]])
  ) # marker order used by joinme_standata()

  list(
    subject = if (active_domains[["subject"]]) {
      build_domain(
        formulae$subject,
        data_event,
        data_long,
        id_variable,
        subject_levels,
        "subject"
      )
    } else {
      list(
        formulae = rep(list(~1), n_classes),
        formula_records = vector("list", n_classes),
        class_specific = FALSE,
        matrix = matrix(0, length(subject_levels), 0L),
        columns = character(0),
        coefficient_class = integer(0),
        coefficient_covariate = character(0),
        class_term_start = rep.int(1L, n_classes),
        class_term_count = integer(n_classes),
        units = as.character(subject_levels)
      )
    },
    marker = if (active_domains[["marker"]]) {
      build_domain(
        formulae$marker,
        data_long,
        NULL,
        marker_variable,
        marker_levels,
        "marker"
      )
    } else {
      list(
        formulae = rep(list(~1), n_classes),
        formula_records = vector("list", n_classes),
        class_specific = FALSE,
        matrix = matrix(0, length(marker_levels), 0L),
        columns = character(0),
        coefficient_class = integer(0),
        coefficient_covariate = character(0),
        class_term_start = rep.int(1L, n_classes),
        class_term_count = integer(n_classes),
        units = as.character(marker_levels)
      )
    },
    class_specific = any(vapply(
      formulae,
      function(specification) isTRUE(specification$class_specific),
      logical(1)
    ))
  )
}

#' Reconstruct one class-regression row for prediction
#'
#' @param design_record Fitted formula, terms, factor levels, contrasts and
#'   column names for one allocation domain.
#' @param primary_data Preferred new-data source.
#' @param fallback_data Optional secondary new-data source.
#' @param unit_variable Column identifying the allocation unit.
#' @param unit_value Identifier of the unit being predicted.
#' @param domain Plain-language domain name used in validation messages.
#'
#' @return A one-row numeric matrix in the fitted design-column order.
#' @keywords internal
#' @noRd
.mixture_prediction_class_design <- function(
  design_record,
  primary_data,
  fallback_data,
  unit_variable,
  unit_value,
  domain
) {
  fitted_columns <- design_record$columns %||% character(0)
  if (length(fitted_columns) == 0L) {
    return(matrix(0, nrow = 1L, ncol = 0L))
  }

  formula_records <- design_record$formula_records # fitted record for every class
  requested_variables <- unique(unlist(
    lapply(formula_records, function(record) {
      all.vars(record$formula)
    }),
    use.names = FALSE
  )) # underlying new-data columns required by any class formula
  unit_frame <- data.frame(
    .class_unit = unit_value,
    stringsAsFactors = FALSE
  ) # single allocation-unit record used to evaluate the fitted formula

  for (variable in requested_variables) {
    source_data <- if (variable %in% names(primary_data)) {
      primary_data
    } else if (
      !is.null(fallback_data) &&
        variable %in% names(fallback_data)
    ) {
      fallback_data
    } else {
      cli::cli_abort(
        "Class covariate {.val {variable}} is unavailable for the new {domain}."
      )
    } # new-data table containing this fitted class covariate

    matching_rows <- as.character(source_data[[unit_variable]]) ==
      as.character(unit_value) # rows belonging to the requested allocation unit
    observed_values <- unique(source_data[[variable]][matching_rows])
    observed_values <- observed_values[!is.na(observed_values)]
    if (length(observed_values) != 1L) {
      cli::cli_abort(c(
        x = "Class covariate {.val {variable}} is not constant for the new {domain}.",
        i = "Supply exactly one non-missing value for each class-regression covariate."
      ))
    }
    unit_frame[[variable]] <- if (is.factor(observed_values)) {
      as.character(observed_values[[1L]])
    } else {
      observed_values[[1L]]
    } # checked covariate value for the new allocation unit
  }

  prediction_designs <- lapply(formula_records, function(record) {
    prediction_frame <- stats::model.frame(
      record$terms,
      data = unit_frame,
      xlev = record$factor_levels %||% list(),
      na.action = stats::na.pass
    ) # formula-ready record using fitted factor levels
    prediction_design <- stats::model.matrix(
      record$terms,
      data = prediction_frame,
      contrasts.arg = record$contrasts %||% NULL
    ) # reconstructed design for this class
    prediction_design <- prediction_design[
      ,
      setdiff(colnames(prediction_design), "(Intercept)"),
      drop = FALSE
    ]
    missing_columns <- setdiff(
      record$columns,
      colnames(prediction_design)
    ) # fitted class columns unavailable in the new record
    if (length(missing_columns) > 0L) {
      cli::cli_abort(
        "The new {domain} cannot reproduce class-design column{?s}: {paste(missing_columns, collapse = ', ')}."
      )
    }
    prediction_design[, record$columns, drop = FALSE]
  }) # class-specific new-data designs in fitted order
  unname(do.call(cbind, prediction_designs))
}

#' Evaluate covariate-dependent latent-class probabilities
#'
#' @param baseline_probability Draw-by-class matrix of baseline probabilities.
#' @param class_coefficient Draw-by-coefficient matrix.
#' @param class_design Unit-by-covariate design matrix.
#' @param class_term_start First coefficient position for every class.
#' @param class_term_count Number of coefficients used by every class.
#'
#' @return Draw-by-unit-by-class array of multinomial-logit probabilities.
#' @keywords internal
#' @noRd
.mixture_class_probability <- function(
  baseline_probability,
  class_coefficient,
  class_design,
  class_term_start,
  class_term_count
) {
  number_draws <- nrow(
    baseline_probability
  ) # retained posterior draws evaluated by dynamic prediction
  n_classes <- ncol(
    baseline_probability
  ) # common number of latent progress classes
  number_units <- nrow(
    class_design
  ) # allocation units whose probabilities are required
  probability <- array(
    0,
    dim = c(number_draws, number_units, n_classes)
  ) # class probabilities to be filled draw by draw and unit by unit

  for (draw in seq_len(number_draws)) {
    baseline_log_probability <- log(
      pmax(baseline_probability[draw, ], .Machine$double.xmin)
    ) # numerically safe fitted baseline logits for this posterior draw
    for (unit in seq_len(number_units)) {
      class_logit <- baseline_log_probability # logits before covariate adjustment
      for (group in seq_len(n_classes)) {
        term_count <- class_term_count[[group]]
        if (term_count > 0L) {
          term_positions <- seq.int(
            class_term_start[[group]],
            length.out = term_count
          ) # concatenated coefficients belonging to this class
          class_logit[[group]] <- class_logit[[group]] +
            sum(
              class_design[unit, term_positions] *
                class_coefficient[draw, term_positions]
            )
        }
      }
      largest_logit <- max(
        class_logit
      ) # centring constant preventing overflow in exponentiation
      unnormalised_probability <- exp(
        class_logit - largest_logit
      ) # stable positive class weights
      probability[draw, unit, ] <- unnormalised_probability /
        sum(unnormalised_probability)
    }
  }
  probability
}

#' Canonicalise latent-progress level names
#'
#' @param class_type User-supplied class types.
#'
#' @return Unique validated names in their original order.
#' @keywords internal
#' @noRd
.canonical_class_types <- function(class_type) {
  supplied_types <- tolower(trimws(
    as.character(class_type %||% character(0))
  )) # normalised public values before checking the deliberately small vocabulary
  accepted_types <- c(
    "subject",
    "marker",
    "corr",
    "vcov"
  ) # terminology shared with the ordinary JoiNMe model
  unknown_types <- setdiff(supplied_types, accepted_types)
  if (length(unknown_types) > 0L) {
    cli::cli_abort(c(
      x = "Unknown {.arg class_type} value{?s}: {paste(unknown_types, collapse = ', ')}.",
      i = "Use only subject, marker, corr, or vcov."
    ))
  }
  validated_types <- unique(supplied_types)
  if (all(c("corr", "vcov") %in% validated_types)) {
    cli::cli_abort(c(
      x = "{.arg class_type} cannot contain both corr and vcov.",
      i = "They are alternative representations of the same marker-by-subject covariance block."
    ))
  }
  validated_types
}

#' Locate correlation coordinates in the packed covariance regression
#'
#' @param q_dimension Number of marker-by-subject random-effect terms.
#'
#' @return One-based positions of the strictly lower-triangular coordinates in
#'   the packed lower triangle used by Stan.
#' @keywords internal
#' @noRd
.mixture_correlation_source_indices <- function(q_dimension) {
  q_dimension <- as.integer(q_dimension) # covariance-matrix order
  packed_position <- 0L # current lower-triangular position
  correlation_positions <- integer(0) # off-diagonal positions returned to Stan
  for (row in seq_len(q_dimension)) {
    for (column in seq_len(row)) {
      packed_position <- packed_position + 1L
      if (column < row) {
        correlation_positions <- c(
          correlation_positions,
          packed_position
        )
      }
    }
  }
  correlation_positions
}

#' Resolve selected coordinates for one class-specific block
#'
#' @param requested User-supplied coordinate indices, or `NULL`.
#' @param available_dimension Dimension of the fitted random-effect block.
#' @param level Canonical block name, used in messages.
#' @param include_all_by_default Whether the default includes every coordinate.
#'
#' @return Integer vector of selected coordinates.
#' @keywords internal
#' @noRd
.mixture_coordinate_indices <- function(
  requested,
  available_dimension,
  level,
  include_all_by_default = FALSE
) {
  available_dimension <- as.integer(available_dimension)
  if (available_dimension < 1L) {
    return(integer(0))
  }
  if (is.null(requested)) {
    if (isTRUE(include_all_by_default)) {
      return(seq_len(available_dimension))
    }
    return(seq_len(min(2L, available_dimension)))
  }

  if (
    !is.numeric(requested) ||
      any(!is.finite(requested)) ||
      any(requested != floor(requested))
  ) {
    cli::cli_abort(c(
      x = "Invalid {.arg class_dimensions} for level {.val {level}}.",
      i = "Choose distinct integer indices between 1 and {available_dimension}."
    ))
  }
  selected <- as.integer(requested)
  if (
    length(selected) == 0L ||
      anyNA(selected) ||
      any(selected < 1L) ||
      any(selected > available_dimension) ||
      anyDuplicated(selected)
  ) {
    cli::cli_abort(c(
      x = "Invalid {.arg class_dimensions} for level {.val {level}}.",
      i = "Choose distinct integer indices between 1 and {available_dimension}."
    ))
  }
  selected
}

#' Find the packed random-intercept coordinate used to identify class labels
#'
#' @param stan_data Prepared ordinary model data and random-effect labels.
#' @param coordinate_indices Selected source coordinates by random-effect block.
#' @param starts First packed coordinate for every block.
#'
#' @return Zero when no selected intercept exists, otherwise its one-based
#'   position in the common mixture vector.
#' @keywords internal
#' @noRd
.mixture_intercept_order_coordinate <- function(
  stan_data,
  coordinate_indices,
  starts
) {
  is_intercept <- function(label) {
    grepl(
      "(^|[^[:alnum:]])intercept([^[:alnum:]]|$)",
      tolower(as.character(label))
    )
  } # whether a fitted random-effect label denotes an intercept
  block_intercepts <- list(
    subject = which(is_intercept(stan_data$zid_cols %||% character(0))),
    marker = which(is_intercept(stan_data$zmk_cols %||% character(0))),
    corr = integer(0),
    vcov = integer(0)
  ) # source coordinates eligible for intercept ordering in each block

  marker_id_labels <- stan_data$zidm_cols %||% character(0)
  marker_id_dimension <- as.integer(stan_data$Q_idm %||% 0L)
  if (marker_id_dimension > 0L && length(marker_id_labels) >= marker_id_dimension) {
    intercept_effects <- which(is_intercept(marker_id_labels))
    if (as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L) {
      block_intercepts$vcov <- intercept_effects
    } else {
      covariance_position <- 0L
      diagonal_positions <- integer(marker_id_dimension)
      for (row in seq_len(marker_id_dimension)) {
        for (column in seq_len(row)) {
          covariance_position <- covariance_position + 1L
          if (row == column) {
            diagonal_positions[[row]] <- covariance_position
          }
        }
      }
      block_intercepts$vcov <- diagonal_positions[intercept_effects]
    }
  }

  for (level in names(coordinate_indices)) {
    selected <- coordinate_indices[[level]]
    intercept_position <- match(
      block_intercepts[[level]],
      selected,
      nomatch = 0L
    ) # local selected position of each eligible intercept
    intercept_position <- intercept_position[intercept_position > 0L]
    if (length(intercept_position) > 0L) {
      return(as.integer(
        starts[[level]] + intercept_position[[1L]] - 1L
      ))
    }
  }
  0L
}

#' Build Stan data for the finite random-effect mixture
#'
#' @description
#' Converts the high-level mixture request into four non-overlapping coordinate
#' segments.  Stan uses these starts and indices to replace the ordinary
#' standardised prior by a finite mixture through an exact log-density ratio.
#'
#' @param stan_data Complete ordinary JoiNMe standata.
#' @param mixture A specification assembled by [joinme_mix()], or `NULL`.
#'
#' @return A named list containing numeric Stan fields and a character metadata
#'   record used only by R reporting methods.
#' @keywords internal
#' @noRd
.build_mixture_standata <- function(stan_data, mixture = NULL) {
  number_subjects <- as.integer(
    stan_data$n_id %||% 0L
  ) # allocation-unit count for the subject class design
  number_markers <- as.integer(
    stan_data$D %||% 0L
  ) # allocation-unit count for the marker class design
  # Ordinary JoiNMe fits use one empty component.  `simplex[1]` has no free
  # parameter, and the zero-column location/scale matrices likewise add no
  # posterior dimension.
  inactive <- is.null(mixture)
  if (inactive) {
    return(list(
      use_mixture = 0L,
      n_classes = 1L,
      K_mix = 0L,
      mix_subject = 0L,
      mix_marker = 0L,
      mix_covariance = 0L,
      mix_dim_subject = 0L,
      mix_idx_subject = integer(0),
      mix_start_subject = 0L,
      mix_dim_marker = 0L,
      mix_idx_marker = integer(0),
      mix_start_marker = 0L,
      mix_dim_covariance = 0L,
      mix_idx_covariance = integer(0),
      mix_start_covariance = 0L,
      mix_ordering = 0L,
      mix_ordered_location_coordinate = 0L,
      mix_probability_prior = 1,
      P_class_subject = 0L,
      X_class_subject = matrix(0, nrow = number_subjects, ncol = 0L),
      class_term_start_subject = 1L,
      class_term_count_subject = 0L,
      P_class_marker = 0L,
      X_class_marker = matrix(0, nrow = number_markers, ncol = 0L),
      class_term_start_marker = 1L,
      class_term_count_marker = 0L,
      prior_class_regression_family = 2L,
      prior_class_regression_mu = numeric(0),
      prior_class_regression_scale = numeric(0),
      prior_class_regression_df = 1,
      prior_class_regression_global_df = 1,
      prior_class_regression_global_scale = 1,
      prior_class_regression_slab_df = 4,
      prior_class_regression_slab_scale = 2,
      mixture = NULL
    ))
  }

  n_classes_requested <- mixture$n_classes
  if (
    !is.numeric(n_classes_requested) ||
      length(n_classes_requested) != 1L ||
      is.na(n_classes_requested) ||
      !is.finite(n_classes_requested) ||
      n_classes_requested < 2 ||
      n_classes_requested != floor(n_classes_requested)
  ) {
    cli::cli_abort("{.arg n_classes} must be a single integer of at least two.")
  }
  n_classes <- as.integer(n_classes_requested)

  selected_types <- .canonical_class_types(mixture$class_type)
  if (length(selected_types) == 0L) {
    cli::cli_abort("{.arg class_type} must select at least one class type.")
  }

  # The covariance-regression block follows the same lower-triangular indexing
  # as Stan's `fit_cov_index.stan`.
  covariance_dimension <- if (
    as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
  ) {
    as.integer(stan_data$Q_idm %||% 0L)
  } else {
    q_dimension <- as.integer(stan_data$Q_idm %||% 0L)
    q_dimension * (q_dimension + 1L) %/% 2L
  }
  correlation_source_indices <- if (
    as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
  ) {
    integer(0)
  } else {
    .mixture_correlation_source_indices(
      as.integer(stan_data$Q_idm %||% 0L)
    )
  } # source positions of correlation rather than scale coordinates
  available_dimensions <- c(
    subject = as.integer(stan_data$R_id %||% 0L),
    marker = as.integer(stan_data$R_mk %||% 0L),
    corr = length(correlation_source_indices),
    vcov = covariance_dimension
  )
  unavailable <- selected_types[available_dimensions[selected_types] < 1L]
  if (length(unavailable) > 0L) {
    reasons <- vapply(unavailable, function(class_type) {
      switch(
        class_type,
        marker = "the longitudinal formula has no marker random-effect coefficients",
        corr = "the formula has no off-diagonal marker-by-subject correlation regression",
        vcov = "the formula has no marker-by-subject covariance regression",
        subject = "the longitudinal formula has no subject random-effect coefficients"
      )
    }, character(1))
    cli::cli_abort(c(
      x = "The requested class type{?s} {?is/are} unavailable: {paste(unavailable, collapse = ', ')}.",
      i = paste(reasons, collapse = "; ")
    ))
  }

  requested_dimensions <- mixture$class_dimensions
  if (!is.null(requested_dimensions) && !is.list(requested_dimensions)) {
    if (length(selected_types) != 1L) {
      cli::cli_abort(
        "{.arg class_dimensions} must be a named list when more than one class type is selected."
      )
    }
    requested_dimensions <- stats::setNames(
      list(requested_dimensions),
      selected_types
    )
  }
  if (is.list(requested_dimensions) && length(requested_dimensions) > 0L) {
    if (is.null(names(requested_dimensions)) || any(!nzchar(names(requested_dimensions)))) {
      cli::cli_abort("{.arg class_dimensions} must be a named list.")
    }
    canonical_dimension_names <- .canonical_class_types(names(requested_dimensions))
    if (anyDuplicated(canonical_dimension_names)) {
      cli::cli_abort("{.arg class_dimensions} supplies the same type more than once.")
    }
    names(requested_dimensions) <- canonical_dimension_names
    extra_dimensions <- setdiff(names(requested_dimensions), selected_types)
    if (length(extra_dimensions) > 0L) {
      cli::cli_abort(
        "{.arg class_dimensions} includes unselected type{?s}: {paste(extra_dimensions, collapse = ', ')}."
      )
    }
  }

  coordinate_indices <- stats::setNames(
    lapply(selected_types, function(class_type) {
      selected_logical_indices <- .mixture_coordinate_indices(
        requested = requested_dimensions[[class_type]] %||% NULL,
        available_dimension = available_dimensions[[class_type]],
        level = class_type
      )
      if (identical(class_type, "corr")) {
        correlation_source_indices[selected_logical_indices]
      } else {
        selected_logical_indices
      }
    }),
    selected_types
  )
  for (level in setdiff(names(available_dimensions), selected_types)) {
    coordinate_indices[[level]] <- integer(0)
  }
  coordinate_indices <- coordinate_indices[names(available_dimensions)]

  dimensions <- vapply(coordinate_indices, length, integer(1))
  starts <- integer(length(dimensions))
  names(starts) <- names(dimensions)
  next_start <- 1L
  for (level in names(dimensions)) {
    if (dimensions[[level]] > 0L) {
      starts[[level]] <- next_start
      next_start <- next_start + dimensions[[level]]
    } else {
      starts[[level]] <- 0L
    }
  }
  total_dimension <- next_start - 1L
  ordering_name <- .resolve_class_ordering(
    class_ordering = mixture$class_ordering,
    class_specific_formulae = isTRUE(mixture$class_design$class_specific)
  ) # validated class-label identification rule
  ordering_code <- match(
    ordering_name,
    c("none", "intercept", "probability")
  ) - 1L # compact Stan code: zero none, one location, two probability
  ordered_location_coordinate <- if (identical(ordering_name, "intercept")) {
    .mixture_intercept_order_coordinate(
      stan_data = stan_data,
      coordinate_indices = coordinate_indices,
      starts = starts
    )
  } else {
    0L
  } # sole packed random-intercept coordinate receiving a strict location order

  probability_prior <- as.numeric(
    mixture$class_probability_concentration %||% 1
  ) # Dirichlet concentration supplied through jm_prior()
  if (length(probability_prior) == 1L) {
    probability_prior <- rep(probability_prior, n_classes)
  }
  if (
    length(probability_prior) != n_classes ||
      any(!is.finite(probability_prior)) ||
      any(probability_prior <= 0)
  ) {
    cli::cli_abort(
      "{.arg priors$class_probability} must be positive and have length one or {.arg n_classes}."
    )
  }
  class_design <- mixture$class_design %||% list(
    subject = list(
      formula = ~1,
      terms = stats::terms(~1),
      factor_levels = list(),
      contrasts = NULL,
      matrix = matrix(0, nrow = number_subjects, ncol = 0L),
      columns = character(0),
      coefficient_class = integer(0),
      coefficient_covariate = character(0),
      class_term_start = rep.int(1L, n_classes),
      class_term_count = integer(n_classes),
      units = as.character(seq_len(number_subjects))
    ),
    marker = list(
      formula = ~1,
      terms = stats::terms(~1),
      factor_levels = list(),
      contrasts = NULL,
      matrix = matrix(0, nrow = number_markers, ncol = 0L),
      columns = character(0),
      coefficient_class = integer(0),
      coefficient_covariate = character(0),
      class_term_start = rep.int(1L, n_classes),
      class_term_count = integer(n_classes),
      units = as.character(seq_len(number_markers))
    )
  ) # checked allocation-unit designs, with an intercept-only compatibility default
  X_class_subject <- class_design$subject$matrix
  X_class_marker <- class_design$marker$matrix
  if (nrow(X_class_subject) != number_subjects) {
    cli::cli_abort(
      "Subject class design does not align with the fitted subject order."
    )
  }
  if (nrow(X_class_marker) != number_markers) {
    cli::cli_abort(
      "Marker class design does not align with the fitted marker order."
    )
  }

  # Materialise the class-regression prior only after both allocation-domain
  # designs have established the exact concatenated coefficient dimension.
  # This delayed expansion permits a concise scalar prior whilst making a
  # coefficient-specific vector unambiguous and safe to validate.
  class_regression_declaration <- mixture$class_regression_prior
  class_regression_declaration <- .normalise_joinme_prior_component(
    class_regression_declaration,
    name = "class_regression",
    default = prior_normal()
  ) # validated family declaration shared by subject and marker class regressions
  class_regression_prior_data <- .materialise_joinme_prior_data(
    priors = list(
      class_regression = .encode_joinme_prior(class_regression_declaration)
    ),
    dimensions = c(
      class_regression = ncol(X_class_subject) + ncol(X_class_marker)
    )
  ) # Stan fields ordered as all subject-domain then all marker-domain coefficients

  include_survival <- isTRUE(mixture$include_survival)

  list(
    use_mixture = 1L,
    n_classes = n_classes,
    K_mix = as.integer(total_dimension),
    mix_subject = as.integer("subject" %in% selected_types),
    mix_marker = as.integer("marker" %in% selected_types),
    mix_covariance = as.integer(any(c("corr", "vcov") %in% selected_types)),
    mix_dim_subject = as.integer(dimensions[["subject"]]),
    mix_idx_subject = as.integer(coordinate_indices$subject),
    mix_start_subject = as.integer(starts[["subject"]]),
    mix_dim_marker = as.integer(dimensions[["marker"]]),
    mix_idx_marker = as.integer(coordinate_indices$marker),
    mix_start_marker = as.integer(starts[["marker"]]),
    mix_dim_covariance = as.integer(
      sum(dimensions[c("corr", "vcov")])
    ),
    mix_idx_covariance = as.integer(c(
      coordinate_indices$corr,
      coordinate_indices$vcov
    )),
    mix_start_covariance = as.integer(max(
      starts[c("corr", "vcov")]
    )),
    mix_ordering = as.integer(ordering_code),
    mix_ordered_location_coordinate =
      as.integer(ordered_location_coordinate),
    mix_probability_prior = probability_prior,
    P_class_subject = as.integer(ncol(X_class_subject)),
    X_class_subject = unname(X_class_subject),
    class_term_start_subject =
      as.integer(class_design$subject$class_term_start),
    class_term_count_subject =
      as.integer(class_design$subject$class_term_count),
    P_class_marker = as.integer(ncol(X_class_marker)),
    X_class_marker = unname(X_class_marker),
    class_term_start_marker =
      as.integer(class_design$marker$class_term_start),
    class_term_count_marker =
      as.integer(class_design$marker$class_term_count),
    prior_class_regression_family = class_regression_prior_data$prior_class_regression_family,
    prior_class_regression_mu = class_regression_prior_data$prior_class_regression_mu,
    prior_class_regression_scale = class_regression_prior_data$prior_class_regression_scale,
    prior_class_regression_df = class_regression_prior_data$prior_class_regression_df,
    prior_class_regression_global_df = class_regression_prior_data$prior_class_regression_global_df,
    prior_class_regression_global_scale = class_regression_prior_data$prior_class_regression_global_scale,
    prior_class_regression_slab_df = class_regression_prior_data$prior_class_regression_slab_df,
    prior_class_regression_slab_scale = class_regression_prior_data$prior_class_regression_slab_scale,
    mixture = list(
      n_classes = n_classes,
      class_type = selected_types,
      dimensions = coordinate_indices,
      starts = starts,
      total_dimension = total_dimension,
      ordering = ordering_name,
      ordered_location_coordinate = ordered_location_coordinate,
      class_probability = probability_prior,
      class_regression_prior = class_regression_declaration,
      class_design = class_design,
      distribution = c(
        "student_t_6",
        "laplace",
        "normal"
      )[as.integer(stan_data$shrinkage %||% 0L) + 1L],
      include_survival = include_survival,
      allocation_domains = list(
        subject = intersect(selected_types, c("subject", "corr", "vcov")),
        marker = intersect(selected_types, "marker")
      )
    )
  )
}
