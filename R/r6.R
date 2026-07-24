#' @name JoiNMe_r6_classes
#' @aliases JoiNMeFit JoiNMeFit JoiNMeDynPred SummaryJoiNMeFit SummaryJoiNMeDynPred
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
#' @export JoiNMeFit JoiNMeFit JoiNMeDynPred SummaryJoiNMeFit SummaryJoiNMeDynPred
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
    print = function(...) {
      .cli_summary_heading("Joint mixed effects model summary", level = 1L)
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

    initialize = function(
      formulaLong,
      formulaEvent,
      formulaVCov,
      parent_call,
      dataLong,
      dataEvent,
      transforms,
      draws,
      threads_per_chain,
      stan_data,
      stan_mod,
      stan_engine,
      stan_args
    ){
      self$formulaLong <- formulaLong
      self$formulaEvent <- formulaEvent
      self$formulaVCov <- formulaVCov
      self$parent_call <- parent_call
      self$dataLong <- dataLong
      self$dataEvent <- dataEvent
      self$transforms <- transforms
      self$draws <- draws
      self$threads_per_chain <- threads_per_chain
      self$stan_data <- stan_data
      self$stan_mod <- stan_mod
      self$stan_engine <- stan_engine
      self$stan_args <- stan_args
    },
    print = function(...) {
      .cli_summary_heading("JoiNMe Stan Data", level = 1L)
      cli::cli_bullets(paste("Stan engine: ", self$stan_engine, "\n", sep = ""))
      cli::cli_bullets("Sample args: ")
      cli::cli_ul(items = 
        if (self$stan_engine == "cmdstanr") {
          self$stan_args
        } else {
          self$stan_args[names(self$stan_args) != "object"]
        }
      )
      invisible(self)
    },
    sample = function(){
      sd <- self$stan_data
      if (self$stan_engine == "cmdstanr") {
        fit <- do.call(self$stan_mod$sample, self$stan_args)
        fit <- .import_cmdstanr_fit(fit)
      } else {
        # Rstan returns warnings about NA in Rhat which is untrue. Can easily work out with diagnosis tho so let's suppress it for now.
        fit <- suppressWarnings(do.call(rstan::sampling, self$stan_args))
      }
      
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
        transform_spec = self$transforms,
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
        dims = c(n_id = sd$n_id, N = sd$N, D = sd$D, P = sd$P, R_id = sd$R_id, R_mk = sd$R_mk, Q_idm = sd$Q_idm),
        draws_default = self$draws,
        threads_per_chain = self$threads_per_chain,
        # tmax_internal = sd$tmax,
        time_indices = list(
          idx_time_beta = sd$idx_time_beta,
          idx_time_uid = sd$idx_time_uid,
          idx_time_vmk = sd$idx_time_vmk,
          idx_time_idm = sd$idx_time_idm
        ),
        marker_weights = sd$marker_weights,
        fixed_marker_weights = sd$fixed_marker_weights,
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

      fit_obj <- JoiNMeFit$new(
        fit = fit,
        stan_data = sd,
        formulaLong = self$formulaLong,
        formulaEvent = self$formulaEvent,
        formulaVCov = self$formulaVCov,
        config = cfg,
        call = self$parent_call,
        tmax = sd$tmax,
        dataLong = self$dataLong,
        dataEvent = self$dataEvent
      )

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
  
