
# ---- corr -----------------------------------------------------------------

#' Extract posterior covariance summaries
#'
#' @description
#' Returns posterior covariance summaries for fitted JoiNMe models. Supports
#' selecting specific blocks through `what` (for example, `what = "id"`).
#'
#' @param object A JoiNMe fit object.
#' @param what Optional covariance block selector (`"id"`, `"marker"`).
#'   When `NULL` (default), returns a nested list for
#'   all covariance components in `formulaLong` and `formulaDist`.
#' @param draws Number of draws to use for summaries.
#' @param ... Unused.
#'
#' @return If `what` is supplied, a data frame for the requested covariance
#'   block. Otherwise, a nested list with `formulaLong` and `formulaDist`
#'   covariance summaries.
#' @method vcov JoiNMeFit
#' @export
vcov.JoiNMeFit <- function(object, what = NULL, draws = NULL, ...) {
 
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!is.null(what)) {
    what <- match.arg(what, choices = c("id", "marker"))
  }

  summarize_cov_matrix <- function(tau_prefix, Lcorr_prefix, dim, label, diagonal_only = FALSE) {
    if (dim <= 0) {
      return(data.frame(
        block = character(0), row = integer(0), col = integer(0),
        Estimate = numeric(0), Est.Error = numeric(0), Q2.5 = numeric(0), Q97.5 = numeric(0),
        Rhat = numeric(0), ess_bulk = numeric(0), ess_tail = numeric(0)
      ))
    }
    tau_vars <- paste0(tau_prefix, "[", seq_len(dim), "]")
    L_vars <- as.vector(outer(seq_len(dim), seq_len(dim),
                              function(r, c) paste0(Lcorr_prefix, "[", r, ",", c, "]")))
    vars <- c(tau_vars, L_vars)

    dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = 1)
    tau_draws <- dmat[, tau_vars, drop = FALSE]
    L_draws <- dmat[, L_vars, drop = FALSE]
    nD <- nrow(tau_draws)

    Sigma_list <- vector("list", nD)
    for (i in seq_len(nD)) {
      tau <- as.numeric(tau_draws[i, ])
      if (isTRUE(diagonal_only)) {
        Corr <- diag(dim)
      } else {
        L <- matrix(as.numeric(L_draws[i, ]), nrow = dim, ncol = dim, byrow = FALSE)
        L[upper.tri(L)] <- 0
        Corr <- L %*% t(L)
      }
      Sigma_list[[i]] <- diag(tau, dim) %*% Corr %*% diag(tau, dim)
    }

    out <- list()
    for (r in seq_len(dim)) for (c in seq_len(dim)) {
      vals <- vapply(Sigma_list, function(S) S[r, c], numeric(1))
      ss <- .summarize_draw_col(vals)
      rhat <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(vals)), error = function(e) NA_real_))
      ess_bulk <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(vals)), error = function(e) NA_real_))
      ess_tail <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(vals)), error = function(e) NA_real_))
      # Keep the same compact inferential schema used in model summaries.
      # Mean, Median, and SD are aliases of Estimate/Est.Error and are
      # intentionally omitted from covariance reporting tables.
      out[[length(out) + 1]] <- data.frame(
        block = label,
        row = r,
        col = c,
        Estimate = as.numeric(ss[["Estimate"]]),
        Est.Error = as.numeric(ss[["Est.Error"]]),
        Q2.5 = as.numeric(ss[["Q2.5"]]),
        Q97.5 = as.numeric(ss[["Q97.5"]]),
        Rhat = rhat,
        ess_bulk = ess_bulk,
        ess_tail = ess_tail,
        stringsAsFactors = FALSE
      )
    }
    df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
    df$row <- as.integer(df$row)
    df$col <- as.integer(df$col)
    df$Estimate <- as.numeric(df$Estimate)
    df$Est.Error <- as.numeric(df$Est.Error)
    df$Q2.5 <- as.numeric(df$Q2.5)
    df$Q97.5 <- as.numeric(df$Q97.5)
    df$Rhat <- as.numeric(df$Rhat)
    df$ess_bulk <- as.numeric(df$ess_bulk)
    df$ess_tail <- as.numeric(df$ess_tail)
    df
  }

  if (!is.null(what)) {
    if (what == "id") {
      return(summarize_cov_matrix(
        "tau_u",
        "Lcorr_u",
        sd$R_id,
        "id",
        diagonal_only = as.integer(sd$indep_id_re %||% 0L) == 1L
      ))
    }
    if (what == "marker") {
      if (sd$R_mk <= 0) {
        cli::cli_abort("Marker corr requested but R_mk = 0.")
      }
      return(summarize_cov_matrix(
        "tau_v",
        "Lcorr_v",
        sd$R_mk,
        "marker",
        diagonal_only = as.integer(sd$indep_marker_re %||% 0L) == 1L
      ))
    }
  }

  summarize_dist_corr <- function(param_name) {
    n_re <- as.integer(sd[[paste0("n_re_", param_name)]] %||% 0L)
    if (n_re <= 0L) return(NULL)

    K_vec <- as.integer(sd[[paste0("K_", param_name)]] %||% integer(0))
    if (length(K_vec) != n_re) return(NULL)

    tau_prefix <- paste0("tau_", param_name)
    term_labels <- object$config$dist$dist_re_terms[[param_name]] %||% rep(NA_character_, n_re)

    out_terms <- vector("list", n_re)
    for (j in seq_len(n_re)) {
      tau_vars <- paste0(tau_prefix, "[", j, ",", seq_len(K_vec[j]), "]")
      dmat <- .get_draws_matrix(fit, variables = tau_vars, draws = draws, seed = 1)
      if (ncol(dmat) == 0) next

      rows <- list()
      idx <- 1L
      for (r in seq_len(K_vec[j])) {
        for (c in seq_len(K_vec[j])) {
          if (r == c) {
            vals <- as.numeric(dmat[, tau_vars[r]])^2
          } else {
            vals <- rep(0, nrow(dmat))
          }
          ss <- .summarize_draw_col(vals)
          rhat <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(vals)), error = function(e) NA_real_))
          ess_bulk <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(vals)), error = function(e) NA_real_))
          ess_tail <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(vals)), error = function(e) NA_real_))
          rows[[idx]] <- data.frame(
            block = term_labels[j] %||% paste0("re_term_", j),
            row = r,
            col = c,
            Estimate = as.numeric(ss[["Estimate"]]),
            Est.Error = as.numeric(ss[["Est.Error"]]),
            Q2.5 = as.numeric(ss[["Q2.5"]]),
            Q97.5 = as.numeric(ss[["Q97.5"]]),
            Rhat = rhat,
            ess_bulk = ess_bulk,
            ess_tail = ess_tail,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
      out_terms[[j]] <- do.call(rbind, rows)
    }
    out_terms <- Filter(Negate(is.null), out_terms)
    if (length(out_terms) == 0) return(NULL)
    out_df <- do.call(rbind, out_terms)
    scope_names <- ifelse(grepl("^family=", out_df$block),
                          sub("^family=([^:]+)::.*$", "\\1", out_df$block),
                          "allFamilies")
    split(out_df, scope_names)
  }

  out <- list(
    formulaLong = list(
      id = summarize_cov_matrix(
        "tau_u",
        "Lcorr_u",
        sd$R_id,
        "id",
        diagonal_only = as.integer(sd$indep_id_re %||% 0L) == 1L
      ),
      marker = if (sd$R_mk > 0) summarize_cov_matrix(
        "tau_v",
        "Lcorr_v",
        sd$R_mk,
        "marker",
        diagonal_only = as.integer(sd$indep_marker_re %||% 0L) == 1L
      ) else NULL
    ),
    formulaDist = list(
      sigma = summarize_dist_corr("sigma"),
      nu = summarize_dist_corr("nu"),
      phi = summarize_dist_corr("phi"),
      alpha = summarize_dist_corr("alpha"),
      kappa = summarize_dist_corr("kappa"),
      tau = summarize_dist_corr("tau")
    )
  )
  out$formulaLong <- out$formulaLong[!vapply(out$formulaLong, is.null, logical(1))]
  out$formulaDist <- out$formulaDist[!vapply(out$formulaDist, is.null, logical(1))]
  out
}


#' Extract predicted covariance summaries
#'
#' @description
#' Returns per-subject covariance summaries for marker-by-id random effects from
#' a `JoiNMeDynPred` object. This method is available only when marker covariance
#' depends on id.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param ... Unused.
#'
#' @return A named list containing covariance summary tables.
#' @method vcov JoiNMeDynPred
#' @export
vcov.JoiNMeDynPred <- function(object, ...) {

    if (!isTRUE(object$metadata$marker_corr_depends_on_id)) {
        cli::cli_abort(c(
            x = "Predicted marker-by-id covariance is only available when marker covariance depends on id.",
            i = "Refit with subject-dependent covariance structure in {.arg formulaVCov} and marker-by-id random effects (Q_idm > 0)."
        ))
    }

    sum_obj <- summary(object)
    list(
        formulaLong = list(
            marker_by_id = sum_obj$tables$corr_marker_id
        )
    )
}
