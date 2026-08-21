#' @name JoiNMe_r6_classes
#' @aliases JoiNMeFit JoiNMeMixFit JoiNMeDynPred JoiNMeMixDynPred SummaryJoiNMeFit SummaryJoiNMeDynPred
#' @title JoiNMe R6 Classes
#'
#' @importFrom R6 R6Class
#'
#' @description
#' R6 classes used to hold JoiNMe fit and prediction results with mutable state.
#'
#' `JoiNMeFit` is the canonical exported fit class. `JoiNMeFit` remains as a
#' compatibility alias.
#'
#' **JoiNMeFit**: holds fitted model, Stan data, and formulas.
#' Public fields: fit, stan_data, formulaLong, formulaEvent, formulaVCov, config, call, tmax,
#' dataLong, dataEvent.
#' Methods: initialize, cache_get, cache_set.
#'
#' **JoiNMeDynPred**: holds dynamic predictions and metadata.
#' Public fields: predictions, quantiles, draws, data, metadata, call, tmax, n_samples.
#' Methods: initialize, cache_get, cache_set.
#'
#' **SummaryJoiNMeFit**: holds cached summary tables and diagnostics.
#' Public fields: tables, diagnostics, metadata.
#' Methods: initialize.
#'
#' **SummaryJoiNMeDynPred**: holds summary tables for predictions.
#' Public fields: tables, metadata.
#' Methods: initialize.
#'
#' @keywords internal
#' @export JoiNMeFit JoiNMeMixFit JoiNMeDynPred JoiNMeMixDynPred SummaryJoiNMeFit SummaryJoiNMeDynPred
NULL


#' @noRd
JoiNMeFit <- R6::R6Class(
  classname = "JoiNMeFit",
  public = list(
    fit = NULL,
    stan_data = NULL,
    formulaLong = NULL,
    formulaEvent = NULL,
    formulaVCov = NULL,
    config = NULL,
    call = NULL,
    tmax = NULL,
    dataLong = NULL,
    dataEvent = NULL,
    #' @description Initialize a `JoiNMeFit` object.
    #' @param fit CmdStanR or rstan fit object.
    #' @param stan_data Stan data list used to fit the model.
    #' @param formulaLong Longitudinal formula used to fit the model.
    #' @param formulaEvent Event formula used to fit the model.
    #' @param formulaVCov Covariance-predictor formula used to fit the model.
    #' @param config Configuration list used to fit the model.
    #' @param call Original function call used to fit the model.
    #' @param tmax Maximum time used to fit the model.
    #' @param dataLong Longitudinal data frame used to fit the model.
    #' @param dataEvent Event data frame used to fit the model.
    #' @return A `JoiNMeFit` object
    initialize = function(
      fit,
      stan_data,
      formulaLong,
      formulaEvent,
      formulaVCov,
      config,
      call,
      tmax,
      dataLong,
      dataEvent
    ) {
      # Store core fit objects and metadata for downstream methods
      self$fit <- fit
      self$stan_data <- stan_data
      self$formulaLong <- formulaLong
      self$formulaEvent <- formulaEvent
      self$formulaVCov <- formulaVCov
      self$config <- config
      self$call <- call
      self$tmax <- tmax
      self$dataLong <- dataLong
      self$dataEvent <- dataEvent
      class(self) <- unique(c("JoiNMeFit", class(self)))
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    },
    cache_clear = function() {
      private$cache <- list()
      invisible(self)
    },
    #' @description Print method for `JoiNMeFit` objects.
    #' @return Invisibly returns `x`.
    print = function(...) {
      .cli_summary_heading(
        if (isTRUE(self$config$include_survival %||% TRUE)) {
          "Joint mixed effects model summary"
        } else {
          "Nested longitudinal mixed effects model summary"
        },
        level = 1L
      )
      if (!is.null(self$call)) {
        cat("Call:\n")
        print(self$call)
      }
      family <- self$stan_data$family_names %||%
        .family_code_to_name(self$config$family_long)
      if (!is.null(family)) {
        fml <-
          if (length(unique(family)) == 1) {
            family[1]
          } else if (length(family) > 5) {
            paste(paste(head(family, 5), collapse = ", "), "...")
          } else {
            paste(family, collapse = ", ")
          }
        cat("Family: ", fml, "\n", sep = "")
      }
      cat(
        "Event process: ",
        if (isTRUE(self$config$include_survival %||% TRUE)) {
          "joint longitudinal-survival"
        } else {
          "longitudinal only"
        },
        "\n",
        sep = ""
      )
      cli::cli_bullets("Use {.code summary()} for parameter summaries.\n")
      invisible(self)
    }
  ),
  private = list(
    cache = list()
  ),
  active = list(
    family = function() {
      sapply(self$config$family_long, .family_name)
    }
  )
)

#' Fitted latent-progress mixture
#'
#' @description
#' R6 container returned by [joinme_mix()].  It inherits every field and method
#' from `JoiNMeFit` and adds the checked mixture description used by
#' class-specific summaries, plots, and predictions.
#'
#' @keywords internal
#' @noRd
JoiNMeMixFit <- R6::R6Class(
  classname = "JoiNMeMixFit",
  inherit = JoiNMeFit,
  public = list(
    mixture = NULL,
    #' @description Initialize a `JoiNMeMixFit` object.
    #' @param fit CmdStanR or rstan fit object.
    #' @param stan_data Stan data list used to fit the model.
    #' @param formulaLong Longitudinal formula used to fit the model.
    #' @param formulaEvent Event formula used to fit the model.
    #' @param formulaVCov Covariance-predictor formula used to fit the model.
    #' @param config Configuration list used to fit the model.
    #' @param call Original function call used to fit the model.
    #' @param tmax Maximum time used to fit the model.
    #' @param dataLong Longitudinal data frame used to fit the model.
    #' @param dataEvent Event data frame used to fit the model.
    #' @param mixture Checked mixture description used to fit the model.
    #' @return A `JoiNMeMixFit` object inheriting from `JoiNMeFit`.
    initialize = function(
      fit,
      stan_data,
      formulaLong,
      formulaEvent,
      formulaVCov,
      config,
      call,
      tmax,
      dataLong,
      dataEvent,
      mixture
    ) {
      super$initialize(
        fit = fit,
        stan_data = stan_data,
        formulaLong = formulaLong,
        formulaEvent = formulaEvent,
        formulaVCov = formulaVCov,
        config = config,
        call = call,
        tmax = tmax,
        dataLong = dataLong,
        dataEvent = dataEvent
      )
      self$mixture <- mixture
      class(self) <- unique(c(
        "JoiNMeMixFit",
        "JoiNMeFit",
        class(self)
      ))
    },
    #' @description Print method for `JoiNMeMixFit` objects.
    #' @return Invisibly returns `x`.
    print = function(...) {
      .cli_summary_heading("Latent-progress joint model", level = 1L)
      if (!is.null(self$call)) {
        cat("Call:\n")
        print(self$call)
      }
      mixture <- self$mixture %||% list()
      cat("Classes: ", mixture$n_classes %||% NA_integer_, "\n", sep = "")
      cat(
        "Class types: ",
        paste(mixture$class_type %||% character(0), collapse = ", "),
        "\n",
        sep = ""
      )
      cat(
        "Event process: ",
        if (isTRUE(mixture$include_survival)) {
          "joint longitudinal-survival"
        } else {
          "longitudinal only"
        },
        "\n",
        sep = ""
      )
      cli::cli_bullets(
        "Use {.code summary()} for class probabilities, locations, and model parameters.\n"
      )
      invisible(self)
    }
  )
)

#' @noRd
JoiNMeDynPred <- R6::R6Class(
  classname = "JoiNMeDynPred",
  public = list(
    predictions = NULL,
    quantiles = NULL,
    draws = NULL,
    data = NULL,
    metadata = NULL,
    call = NULL,
    tmax = NULL,
    n_samples = NULL,
    initialize = function(
      predictions,
      quantiles,
      draws,
      data,
      metadata,
      call,
      tmax,
      n_samples
    ) {
      # Store prediction outputs and metadata for plotting/summaries
      self$predictions <- predictions
      self$quantiles <- quantiles
      self$draws <- draws
      self$data <- data
      self$metadata <- metadata
      self$call <- call
      self$tmax <- tmax
      self$n_samples <- n_samples
      class(self) <- unique(c("JoiNMeDynPred", class(self)))
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    },
    cache_clear = function() {
      private$cache <- list()
      invisible(self)
    }
  ),
  private = list(
    cache = list()
  )
)

#' Dynamic predictions from a latent-progress mixture
#'
#' @description
#' Prediction container inheriting the complete `JoiNMeDynPred` interface and
#' retaining the mixture definition of its parent fit.
#'
#' @noRd
JoiNMeMixDynPred <- R6::R6Class(
  classname = "JoiNMeMixDynPred",
  inherit = JoiNMeDynPred,
  public = list(
    mixture = NULL,

    initialize = function(
      predictions,
      quantiles,
      draws,
      data,
      metadata,
      call,
      tmax,
      n_samples,
      mixture
    ) {
      super$initialize(
        predictions = predictions,
        quantiles = quantiles,
        draws = draws,
        data = data,
        metadata = metadata,
        call = call,
        tmax = tmax,
        n_samples = n_samples
      )
      self$mixture <- mixture
      class(self) <- unique(c(
        "JoiNMeMixDynPred",
        "JoiNMeDynPred",
        class(self)
      ))
    }
  )
)

#' @noRd
SummaryJoiNMeFit <- R6::R6Class(
  classname = "summary_JoiNMeFit",
  public = list(
    tables = NULL,
    diagnostics = NULL,
    metadata = NULL,
    initialize = function(tables, diagnostics, metadata) {
      # Store summary tables and diagnostics for printing
      self$tables <- tables
      self$diagnostics <- diagnostics
      self$metadata <- metadata
      class(self) <- unique(c(
        "summary_JoiNMeFit",
        "summary_JoiNMeFit",
        class(self)
      ))
    }
  )
)

#' @noRd
SummaryJoiNMeDynPred <- R6::R6Class(
  classname = "summary_JoiNMeDynPred",
  public = list(
    tables = NULL,
    metadata = NULL,
    initialize = function(tables, metadata) {
      # Store prediction summary tables
      self$tables <- tables
      self$metadata <- metadata
      class(self) <- unique(c(
        # "summary_PredJoiNMeFit",
        "summary_JoiNMeDynPred",
        class(self)
      ))
    }
  )
)

#' Intermediate standata and config for JoiNMe model fitting
#' @keywords internal
#' @noRd
JoiNMeStanData <- R6::R6Class(
  classname = "JoiNMeStanData",
  public = list(
    formulaLong = NULL,
    formulaEvent = NULL,
    formulaVCov = NULL,
    parent_call = NULL,
    dataLong = NULL,
    dataEvent = NULL,
    transforms = NULL,
    draws = NULL,
    threads_per_chain = NULL,
    stan_data = NULL,
    stan_mod = NULL,
    stan_engine = NULL,
    stan_args = NULL,
    result_class = NULL,
    result_metadata = NULL,

    initialize = function(
      recipes
    ) {
      self$formulaLong <- recipes$formulaLong
      self$formulaEvent <- recipes$formulaEvent
      self$formulaVCov <- recipes$formulaVCov
      self$parent_call <- recipes$parent_call
      self$dataLong <- recipes$dataLong
      self$dataEvent <- recipes$dataEvent
      self$transforms <- recipes$transforms
      self$draws <- recipes$draws
      self$threads_per_chain <- recipes$threads_per_chain
      self$stan_data <- recipes$stan_data
      self$stan_mod <- recipes$stan_mod
      self$stan_engine <- recipes$stan_engine
      self$stan_args <- recipes$stan_args
    },
    print = function(...) {
      .cli_summary_heading("JoiNMe Stan Data", level = 1L)
      cli::cli_bullets(paste("Stan engine: ", self$stan_engine, "\n", sep = ""))
      # cli::cli_bullets("Sample args: ")
      # cli::cli_ul(
      #   items = if (self$stan_engine == "cmdstanr") {
      #     self$stan_args
      #   } else {
      #     self$stan_args[names(self$stan_args) != "object"]
      #   }
      # )
      invisible(self)
    },
    make_cfg = function() {
      sd <- self$stan_data
      cfg <- list(
        family_long = sd$family_long,
        family_link = sd$link_long,
        family_link_names = sd$link_names,
        family_inv_link_n_ops = sd$inv_link_n_ops,
        family_inv_link_ops = sd$inv_link_ops,
        family_inv_link_n_const = sd$inv_link_n_const,
        family_inv_link_const = sd$inv_link_const,
        assoc = c(
          cv_total = sd$assoc_cv_total,
          cv_mean = sd$assoc_cv_mean,
          cv_marker = sd$assoc_cv_marker,
          cs_total = sd$assoc_cs_total,
          cs_mean = sd$assoc_cs_mean,
          cs_marker = sd$assoc_cs_marker,
          corr = sd$assoc_corr,
          vcov = sd$assoc_vcov
        ),
        transforms = list(
          tf_mode_cv_tot = sd$tf_mode_cv_tot,
          tf_mode_cs_tot = sd$tf_mode_cs_tot,
          tf_mode_cv_mean = sd$tf_mode_cv_mean,
          tf_mode_cs_mean = sd$tf_mode_cs_mean,
          tf_mode_cv_marker = sd$tf_mode_cv_marker,
          tf_mode_cs_marker = sd$tf_mode_cs_marker,
          tf_mode_corr = sd$tf_mode_corr,
          tf_mode_vcov = sd$tf_mode_vcov
        ),
        # Retain both names for objects written by earlier development
        # versions.  Plotting and class-trajectory reconstruction use the
        # plural field as the canonical resolved specification.
        transform_spec = self$transforms,
        transforms_spec = self$transforms,
        dist = list(
          dist_cols = sd$dist_cols,
          dist_re_terms = sd$dist_re_terms,
          dist_formulas = sd$dist_formulas
        ),
        indep = c(
          id = sd$indep_id_re,
          marker = sd$indep_marker_re,
          idmarker_cov = sd$indep_idmarker_cov
        ),
        allow_marker_crosscorr = sd$allow_marker_crosscorr,
        shrinkage = sd$shrinkage,
        dims = c(
          n_id = sd$n_id,
          N = sd$N,
          D = sd$D,
          P = sd$P,
          R_id = sd$R_id,
          R_mk = sd$R_mk,
          Q_idm = sd$Q_idm
        ),
        draws_default = self$draws,
        threads_per_chain = self$threads_per_chain,
        include_survival = isTRUE(
          as.integer(sd$include_survival %||% 1L) == 1L
        ), # whether survival observations contributed to the fitted likelihood
        marker_weight_offsets = sd$marker_weight_offsets,
        basehaz = sd$basehaz,
        n_knots = sd$basehaz_n_knots,
        basehaz_n_knots = sd$basehaz_n_knots,
        basehaz_knots = sd$basehaz_knots,
        basehaz_degree = sd$basehaz_degree,
        basehaz_formula = sd$basehaz_formula,
        basehaz_col_means = sd$basehaz_col_means,
        K_event = sd$K_event,
        surv_type = sd$surv_type,
        event_censor_types_present = sd$event_censor_types_present,
        vcov_diag_link = sd$vcov_diag_link,
        use_tau_fixed = sd$use_tau_fixed,
        tau_fixed = sd$tau_fixed
      )
      cfg$engine <- self$stan_engine
      cfg
      # cfg$mixture <- self$result_metadata %||% sd$mixture
    },
    sample = function(...) {
      sd <- self$stan_data
      if (length(list(...)) > 0) {
        self$stan_args <- modifyList(self$stan_args, list(...))
      }

      if (self$stan_engine == "cmdstanr") {
        fit <- do.call(self$stan_mod$sample, self$stan_args)
        fit <- .import_cmdstanr_fit(fit)
      } else {
        # Rstan returns warnings about NA in Rhat which is untrue. Can easily work out with diagnosis tho so let's suppress it for now.
        fit <- suppressWarnings(do.call(rstan::sampling, self$stan_args))
      }

      fit_obj <- JoiNMeFit$new(
        fit = fit,
        stan_data = sd,
        formulaLong = self$formulaLong,
        formulaEvent = self$formulaEvent,
        formulaVCov = self$formulaVCov,
        config = self$make_cfg(),
        call = self$parent_call,
        tmax = sd$tmax,
        dataLong = self$dataLong,
        dataEvent = self$dataEvent
      )

      # fit_obj <- if (identical(self$result_class, "mixture")) {
      #   do.call(
      #     JoiNMeMixFit$new,
      #     modifyList(
      #       common_fit_arguments,
      #       list(mixture = self$result_metadata %||% sd$mixture)
      #     )
      #   )
      # } else {
      #   do.call(JoiNMeFit$new, common_fit_arguments)
      # }

      # Store a plotting bundle for association plots remain
      # usable even when cmdstanr CSV outputs are no longer available.
      # fit_obj$config$association_plot_ <- tryCatch(
      #   .build_JoiNMefit_association_plot_(
      #     fit = fit,
      #     stan_data = sd,
      #     config = fit_obj$config,
      #     dataLong = self$dataLong,
      #     seed = self$stan_args$seed %||% defaults$seed
      #   ),
      #   error = function(e) NULL
      # )

      fit_obj
    }
  )
)

JoiNMeMixStanData <- R6::R6Class(
  classname = "JoiNMeMixStanData",
  inherit = JoiNMeStanData,
  public = list(
    mixture = NULL,
    initialize = function(
      recipe,
      more_recipe = list()
    ) {
      super$initialize(recipe)
      self$mixture <- more_recipe$mixture
    },
    print = function() {
      .cli_summary_heading("Latent mixture Stan Data", level = 1L)
      cli::cli_bullets(paste("Stan engine: ", self$stan_engine, "\n", sep = ""))
      # cli::cli_bullets("Sample args: ")
      # cli::cli_ul(
      #   items = if (self$stan_engine == "cmdstanr") {
      #     self$stan_args
      #   } else {
      #     self$stan_args[names(self$stan_args) != "object"]
      #   }
      # )
      cli::cli_bullets(paste(
        "Classes: ",
        self$mixture$n_classes %||% NA_integer_,
        "\n",
        sep = ""
      ))
      cli::cli_bullets(paste(
        "Class types: ",
        paste(self$mixture$class_type %||% character(0), collapse = ", "),
        "\n",
        sep = ""
      ))
      cli::cli_bullets(paste(
        "Event process: ",
        if (isTRUE(self$mixture$include_survival)) {
          "joint longitudinal-survival"
        } else {
          "longitudinal only"
        },
        "\n",
        sep = ""
      ))
      invisible(self)
    },
    make_cfg = function() {
      cfg <- super$make_cfg()
      cfg$mixture <- self$mixture
      cfg
    },
    sample = function() {
      sd <- self$stan_data
      if (self$stan_engine == "cmdstanr") {
        fit <- do.call(self$stan_mod$sample, self$stan_args)
        fit <- .import_cmdstanr_fit(fit)
      } else {
        # Rstan returns warnings about NA in Rhat which is untrue. Can easily work out with diagnosis tho so let's suppress it for now.
        fit <- suppressWarnings(do.call(rstan::sampling, self$stan_args))
      }

      fit_obj <- JoiNMeMixFit$new(
        fit = fit,
        stan_data = sd,
        formulaLong = self$formulaLong,
        formulaEvent = self$formulaEvent,
        formulaVCov = self$formulaVCov,
        config = self$make_cfg(),
        call = self$parent_call,
        tmax = sd$tmax,
        dataLong = self$dataLong,
        dataEvent = self$dataEvent,
        mixture = self$mixture
      )

      fit_obj
    }
  )
)
