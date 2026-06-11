#' Enhanced Plot Method for Dynamic Prediction from Joint Models
#'
#' @importFrom stats na.omit
#' @importFrom stats setNames
#'
#' @description
#' Produces publication-ready ggplot2 visualisations of dynamic predictions from
#' a joint model with flexible customisation of longitudinal trajectories and
#' conditional survival curves.
#'
#' @param x An object of class `JoinMeDynPred` returned by the [predict] method.
#' @param type Character vector indicating what to plot: "cumhaz" (default),
#'   "longitudinal", or "survival". Multiple values are allowed.
#' @param subject Integer/Character vector. Which subject(s) to plot. If NULL,
#'   plots all subjects in separate panels or returns a list.
#' @param marker Optional marker subset for longitudinal plots. Use a specific
#'   marker level to plot only that marker, or `NA` to plot all markers.
#' @param trajectory_type Character. Which trajectory to show:
#'   "id_marker" (default), "marker_pop" (marker-level average),
#'   "overall_pop" (overall mean), or "all_types" for all three.
#' @param scale Optional longitudinal prediction scale to plot. Must match
#'   available scales in the prediction object ("epred", "linpred", "predict").
#'   If NULL, defaults to "epred" when available, otherwise the first
#'   available scale.
#' @param smooth_trajectory Logical. Whether to smooth trajectories with splines
#'   (default TRUE). Only applicable if sufficient observations.
#' @param smooth_method Character. Smoothing method: "loess" or "spline".
#' @param smooth_span Numeric in (0,1). Controls span/flexibility of LOESS smoothing
#'   (default 0.3).
#' @param ci_levels Numeric vector. Credible interval levels to display
#'   (default c(0.5, 0.95) for 50% and 95%).
#' @param ci_type Character. How to display uncertainty: "ribbon" (shaded area),
#'   "line" (quantile lines), or "both".
#' @param observed_first Logical. Whether to visually separate observed history
#'   (shaded background) from prediction region (default TRUE).
#' @param facet_by Character. "marker" to facet by marker type, "none" for single plot.
#' @param facet_scales Character. "fixed" or "free_y" for longitudinal faceting.
#' @param combined Logical. If TRUE and multiple plot types are requested, combine
#'   them into a single layout (patchwork/cowplot).
#'   - For a single subject, returns one combined plot.
#'   - For multiple subjects, returns a named list of combined plots (one per subject).
#'   - If a combiner package is unavailable, falls back to the uncombined
#'     subject/outcome structure for that subject.
#' @param show_data Logical. Whether to overlay observed data points (default TRUE).
#' @param show_observed_line Logical. Whether to connect observed points with lines
#'   (default TRUE).
#' @param observed_style List with elements "color", "shape", "size", "alpha" for
#'   appearance of observed data points.
#' @param prediction_style List with elements "color", "fill", "linewidth", "alpha"
#'   for prediction line and ribbons.
#' @param theme_fn ggplot2 theme function (default `ggplot2::theme_minimal`).
#' @param palette_marker Character vector of colors for different markers, or
#'   function like `ggplot2::scale_color_brewer()`.
#' @param ... Additional arguments (unused).
#'
#' @details
#' **Longitudinal Plots:**
#' - The longitudinal scale is selected via `.arg scale` (or inferred when NULL).
#' - On \"predict\" scale: observed values are shown for times in interval
#'   `[0, time_start]` and predictions are drawn for `[time_start, time_horizon]`.
#' - On \"linpred\" or \"epred\" scale: predictions are shown for `[0, time_horizon]`,
#'   with fitted values providing the pre-time_start segment when available.
#' - Optional shaded region distinguishing observed from prediction periods.
#' - Predicted trajectories with credible bands.
#' - Multiple CI levels with decreasing alpha for visual hierarchy.
#' - Optional smoothing for smoother appearance.
#' - Support for multiple markers with automatic faceting.
#'
#' **Survival / Cumulative Hazard Plots:**
#' - Conditional survival probability S(t | T_cond).
#' - Conditional cumulative hazard H(t | T_cond) (default).
#' - Credible bands with multiple levels.
#'
#' @return
#' If a single subject and outcome requested: a `ggplot` object.
#' If multiple subjects and `combined = TRUE`: a named list where each element is
#'   a per-subject combined plot (patchwork/cowplot) when backend support exists,
#'   otherwise the corresponding uncombined per-subject list.
#' If multiple subjects and `combined = FALSE`: a named list grouped by subject,
#'   each containing requested outcome plots.
#' If multiple outcomes: list with "longitudinal", "survival", and/or "cumhaz" elements,
#'   or a single combined plot if `combined = TRUE`.
#'
#' @export
# File overview:
# - Build longitudinal and survival/cumhaz plots from prediction summaries.
# - Support per-subject panels and combined outputs.
plot.JoinMeDynPred <- function(
    x,
    type = c("cumhaz", "longitudinal", "survival"),
    subject = NULL,
    marker = NA,
    trajectory_type = c("id_marker", "marker_pop", "overall_pop", "all_types"),
    scale = NULL,
    smooth_trajectory = TRUE,
    smooth_method = c("loess", "spline"),
    smooth_span = 0.3,
    ci_levels = c(0.5, 0.95),
    ci_type = c("ribbon", "line", "both"),
    observed_first = TRUE,
    facet_by = c("marker", "none"),
    facet_scales = "free_y",
    combined = TRUE,
    show_data = TRUE,
    show_observed_line = TRUE,
    observed_style = list(
        color = "black",
        shape = 21,
        size = 2,
        alpha = 0.6
    ),
    prediction_style = list(
        color = "steelblue",
        fill = "steelblue",
        linewidth = 0.8,
        alpha = 0.2
    ),
    theme_fn = ggplot2::theme_bw,
    palette_marker = NULL,
    ...
) {
    dots <- list(...)
    is_competing <- FALSE
    if (!is.null(x$data$event)) {
        ev_cols <- names(x$data$event)
        is_competing <- any(c("event_type", "stage_from", "stage_to") %in% ev_cols)
    }
    if ("which" %in% names(dots)) {
        type <- dots$which
    }
    if (missing(type) && !("which" %in% names(dots))) {
        type <- if (isTRUE(is_competing)) "cumhaz" else c("longitudinal", "survival")
    }
    type <- .validate_plot_type(type)
    trajectory_type <- match.arg(trajectory_type)
    smooth_method <- match.arg(smooth_method)
    facet_by <- match.arg(facet_by)
    ci_type <- match.arg(ci_type)

    ci_levels <- .validate_ci_levels_plot(ci_levels, x$metadata$ci_levels %||% NULL)
    marker <- .normalize_plot_marker_filter(marker, available_markers = .available_plot_markers(x))

    available_scales <- .available_longitudinal_scales(x)
    if (!is.null(scale)) {
        scale <- match.arg(scale, choices = available_scales)
    } else if (length(available_scales) > 0) {
        scale <- if ("epred" %in% available_scales) "epred" else available_scales[1]
    }

    # Identify subjects
    # - allow explicit subject list or infer from prediction summaries
    ids <- if (!is.null(subject)) subject else {
        unique(c(
            x$quantiles$longitudinal$id %||% x$predictions$longitudinal$id,
            x$quantiles$survival$id %||% x$predictions$survival$id,
            x$quantiles$cumhaz$id %||% x$predictions$cumhaz$id
        ))
    }
    ids <- na.omit(ids)
    if (length(ids) == 0) {
        cli::cli_warn(c(
            x = "No subjects found for plotting.",
            i = "Provide a valid {.arg subject} or check the prediction object."
        ))
        return(invisible(NULL))
    }

    plots_list <- list()
    
    # Determine trajectory sources based on type
    # - multiple trajectories when trajectory_type == "all_types"
    trajectory_sources <- if (trajectory_type == "all_types") {
        c("id_marker", "marker_pop", "overall_pop")
    } else {
        trajectory_type
    }

    for (id in ids) {
        plots_sub <- list()

        # =====================================================================
        # LONGITUDINAL PLOT
        # =====================================================================
        if ("longitudinal" %in% type) {
            for (traj_src in trajectory_sources) {
                p_long <- .plot_longitudinal_single(
                    x, id, traj_src,
                    marker = marker,
                    scale_long = scale,
                    smooth_trajectory, smooth_method, smooth_span,
                    ci_levels, ci_type, observed_first,
                    facet_by, facet_scales,
                    show_data, show_observed_line,
                    observed_style, prediction_style,
                    theme_fn, palette_marker
                )
                if (length(trajectory_sources) > 1) {
                    plots_sub[[paste0("longitudinal_", traj_src)]] <- p_long
                } else {
                    plots_sub$longitudinal <- p_long
                }
            }
        }

        # =====================================================================
        # SURVIVAL PLOT
        # =====================================================================
        if ("survival" %in% type) {
            p_surv <- .plot_survival_single(
                x, id,
                ci_levels, ci_type,
                prediction_style, theme_fn
            )
            plots_sub$survival <- p_surv
        }

        # =====================================================================
        # CUMHAZ PLOT
        # =====================================================================
        if ("cumhaz" %in% type) {
            p_cumhaz <- .plot_cumhaz_single(
                x, id,
                ci_levels, ci_type,
                prediction_style, theme_fn
            )
            plots_sub$cumhaz <- p_cumhaz
        }

        plots_list[[as.character(id)]] <- plots_sub
    }

    # =========================================================================
    # RETURN
    # =========================================================================
        if (length(plots_list) == 1) {
            # Single subject
            res <- plots_list[[1]]
            if (length(res) == 1) {
                return(res[[1]])
            }
            if (isTRUE(combined)) {
                combined_plot <- .combine_plot_grid(res, ncol = length(type))
                if (inherits(combined_plot, "gg")) return(combined_plot)
            }
            return(res)
        }

        # Multiple subjects
        if (isTRUE(combined)) {
            # Combine requested outcomes within each subject.
            # Output structure for multiple subjects:
            # - preferred: list(id -> combined plot object)
            # - fallback: list(id -> uncombined per-subject plot list)
            out <- lapply(plots_list, function(sub_plots) {
                design <- switch(
                    length(type),
                    "1" = 'A',
                    "2" = "AB",
                    "3" = "AB\nAC"
                )
                .combine_plot_grid(
                    sub_plots,
                    ncol = min(2, length(type)),
                    nrow = min(1, length(type) - 1),
                    design = design,
                    fallback = "input"
                )
            })
            return(out)
        }
        if (length(type) == 1 && type == "longitudinal") {
            return(lapply(plots_list, function(x) x$longitudinal))
        }
        if (length(type) == 1 && type == "survival") {
            return(lapply(plots_list, function(x) x$survival))
        }
        if (length(type) == 1 && type == "cumhaz") {
            return(lapply(plots_list, function(x) x$cumhaz))
        }
        plots_list
}

# ============================================================================
# Helper: Longitudinal Plot for Single Subject
# ============================================================================
.plot_longitudinal_single <- function(
    x, id, trajectory_type = "id_marker",
    marker = NULL,
    scale_long,
    smooth_trajectory, smooth_method, smooth_span,
    ci_levels, ci_type, observed_first,
    facet_by, facet_scales,
    show_data, show_observed_line,
    observed_style, prediction_style,
    theme_fn, palette_marker
) {

    # Extract quantile data based on trajectory type
    # - population vs subject-specific paths
    if (trajectory_type == "marker_pop") {
        quant_df <- x$quantiles$longitudinal_marker_pop
        if (is.null(quant_df) || nrow(quant_df) == 0) {
            cli::cli_warn(c(
                x = "Population-level marker trajectories not available.",
                i = "Check that the prediction output includes marker-population summaries."
            ))
            return(NULL)
        }
        show_data <- FALSE  # Don't show individual observed data for population
    } else if (trajectory_type == "overall_pop") {
        quant_df <- x$quantiles$longitudinal_overall_pop
        if (is.null(quant_df) || nrow(quant_df) == 0) {
            cli::cli_warn(c(
                x = "Overall population trajectory not available.",
                i = "Check that the prediction output includes overall-population summaries."
            ))
            return(NULL)
        }
        show_data <- FALSE  # Don't show individual observed data for population
    } else {
        quant_df <- x$quantiles$longitudinal
    }
    
    if (!is.null(scale_long) && !is.null(quant_df) && "scale" %in% names(quant_df)) {
        quant_df <- quant_df[quant_df$scale == scale_long, , drop = FALSE]
    }

    if (is.null(quant_df)) quant_df <- x$predictions$longitudinal
    if (!is.null(scale_long) && !is.null(quant_df) && "scale" %in% names(quant_df)) {
        quant_df <- quant_df[quant_df$scale == scale_long, , drop = FALSE]
    }
    if (is.null(quant_df)) {
        cli::cli_warn(c(
            x = "No quantile or prediction data available for subject {id} ({trajectory_type}).",
            i = "Ensure the prediction includes the requested trajectory type."
        ))
        return(NULL)
    }

    # Filter by id if not population level
    if (trajectory_type == "id_marker") {
        quant_df <- quant_df[quant_df$id == id, , drop = FALSE]
    } else if (trajectory_type == "marker_pop") {
        # Already filtered - use all rows
    } else if (trajectory_type == "overall_pop") {
        # Already filtered - use all rows
    }
    if (!is.null(marker) && "marker" %in% names(quant_df)) {
        quant_df <- quant_df[as.character(quant_df$marker) %in% marker, , drop = FALSE]
    }
    
    if (nrow(quant_df) == 0) {
        cli::cli_warn(c(
            x = "No data found for subject {id} ({trajectory_type}).",
            i = "Check the prediction output and subject identifiers."
        ))
        return(NULL)
    }
    if (is.null(scale_long)) {
        scale_long <- x$metadata$scale %||% "epred"
    }
    if (!scale_long %in% c("epred", "linpred", "predict")) scale_long <- "epred"

    # Get observed data
    data_long <- x$data$longitudinal
    data_long_id <- data_long[data_long$id == id, , drop = FALSE]

    # Identify variable names from metadata (stored during prediction)
    # Fall back to x$call if not available
    id_var <- x$metadata$id_var %||% eval(x$call$id_var) %||% "id"
    time_var <- x$metadata$time_var %||% eval(x$call$time_var) %||% "time"
    marker_var <- x$metadata$marker_var %||% eval(x$call$marker_var) %||% "marker"
    if (!is.null(marker) && marker_var %in% names(data_long_id)) {
        data_long_id <- data_long_id[as.character(data_long_id[[marker_var]]) %in% marker, , drop = FALSE]
    }
    resp_var <- x$metadata$response_var %||% NA_character_
    
    # Ensure resp_var is actually a character string, not a formula object or invalid name
    if (!is.character(resp_var) || length(resp_var) == 0 || is.na(resp_var)) {
        resp_var <- NA_character_
    } else if (!(resp_var %in% names(data_long_id))) {
        # If resp_var is not a valid column, treat as NA and fall back
        resp_var <- NA_character_
    }
    
    # If response variable still not found, try to extract it from formula
    if (is.na(resp_var)) {
        formula_to_use <- if (!is.null(x$call$formulaLong)) x$call$formulaLong else NA
        if (!identical(formula_to_use, NA)) {
            resp_var <- tryCatch(
                all.vars(formula_to_use)[1],
                error = function(e) NA_character_
            )
            # Validate the extracted variable exists in data
            if (!is.na(resp_var) && !(resp_var %in% names(data_long_id))) {
                resp_var <- NA_character_
            }
        }
    }
    
    # Final fallback: try to guess from data if still missing
    if (is.na(resp_var)) {
        # Check if data has a column that looks like response
        resp_candidates <- setdiff(names(data_long_id), c(id_var, time_var, marker_var))
        # Prioritize "y" as response variable if it exists
        if ("y" %in% resp_candidates) {
            resp_var <- "y"
        } else if (length(resp_candidates) > 0) {
            # Otherwise take the first non-key column
            resp_var <- resp_candidates[1]
        } else {
            cli::cli_abort(c(
                x = "Cannot determine response variable name for longitudinal data.",
                i = "Ensure the prediction object includes metadata$response_var or a valid formulaLong with a response, or that the data contains a recognizable response column."
            ))
        }
    }

    # Conditioning time
    # - separates observed history from prediction horizon
    t_cond <- .conditioning_time_for_id(x, id)

    # Initialise plot data
    df_pred <- quant_df
    df_pred$type <- "prediction"
    df_pred$segment <- "predict"
    df_pred$marker <- as.character(df_pred$marker)

    # Prepare observed or fitted data for the left-of-time_start region
    df_obs <- NULL
    use_fitted <- scale_long %in% c("linpred", "epred")

    if (isTRUE(show_data)) {
        if (use_fitted) {
            fit_df <- x$quantiles$longitudinal_fitted
            if (!is.null(fit_df) && "scale" %in% names(fit_df)) {
                fit_df <- fit_df[fit_df$id == id & fit_df$scale == scale_long, , drop = FALSE]
                if (!is.null(marker) && "marker" %in% names(fit_df)) {
                    fit_df <- fit_df[as.character(fit_df$marker) %in% marker, , drop = FALSE]
                }
                if (nrow(fit_df) > 0) {
                    median_col <- .quantile_name_from_prob(0.5)
                    if (!(median_col %in% names(fit_df))) median_col <- "mean"
                    df_obs <- fit_df[, c("time", "marker", median_col), drop = FALSE]
                    names(df_obs)[names(df_obs) == median_col] <- "value"
                    df_obs$type <- ifelse(df_obs$time <= t_cond, "fitted", "fitted_future")
                    df_obs$marker <- as.character(df_obs$marker)
                }
            }
            if (is.null(df_obs) || nrow(df_obs) == 0) df_obs <- NULL
        } else if (nrow(data_long_id) > 0) {
            # Verify that required columns exist in the data
            required_cols <- c(time_var, marker_var, resp_var)
            missing_cols <- setdiff(required_cols, names(data_long_id))
            if (length(missing_cols) > 0) {
                cli::cli_warn(c(
                    x = "Cannot plot observed data: missing columns: {paste(missing_cols, collapse = ', ')}.",
                    i = "Check the prediction input data and column names."
                ))
            } else {
                df_obs <- data_long_id[, c(time_var, marker_var, resp_var), drop = FALSE]
                names(df_obs) <- c("time", "marker", "value")
                df_obs$type <- ifelse(df_obs$time <= t_cond, "observed", "observed_future")
                df_obs$marker <- as.character(df_obs$marker)
            }
        }
    }

    # Build prediction coverage from 0 to horizon
    t_horiz <- max(quant_df$time, na.rm = TRUE)
    if (!is.finite(t_horiz) && !is.null(df_obs)) {
        t_horiz <- max(df_obs$time, na.rm = TRUE)
    }
    if (!is.finite(t_horiz)) t_horiz <- t_cond

    if (scale_long == "predict") {
        df_pred <- df_pred[df_pred$time >= t_cond & df_pred$time <= t_horiz, , drop = FALSE]
        if (!is.null(df_obs)) df_obs <- df_obs[df_obs$time <= t_cond, , drop = FALSE]
    } else {
        df_pred_future <- df_pred[df_pred$time >= t_cond & df_pred$time <= t_horiz, , drop = FALSE]
        df_pred_future$segment <- "predict"
        df_pred_hist <- NULL
        fit_quant <- x$quantiles$longitudinal_fitted
        if (!is.null(fit_quant)) {
            fit_quant <- fit_quant[fit_quant$id == id & fit_quant$scale == scale_long, , drop = FALSE]
            if (!is.null(marker) && "marker" %in% names(fit_quant)) {
                fit_quant <- fit_quant[as.character(fit_quant$marker) %in% marker, , drop = FALSE]
            }
            if (nrow(fit_quant) > 0) {
                df_pred_hist <- fit_quant
                df_pred_hist$segment <- "fitted"
                df_pred_hist <- df_pred_hist[df_pred_hist$time <= t_cond, , drop = FALSE]
            }
        }
        df_pred <- if (!is.null(df_pred_hist)) {
            tidytable::bind_rows(df_pred_hist, df_pred_future)
        } else {
            df_pred_future
        }
    }

    # Start plot with data and base aesthetics
    if (is.null(df_pred) || nrow(df_pred) == 0) {
        stop("No valid data to plot for subject ", id)
    }
    
    # Initialise plot - don't set y in base aes since different layers use different y columns
    p <- ggplot2::ggplot(data = df_pred, ggplot2::aes(x = time, fill = marker, color = marker))

    # Add observed region background (if requested)
    if (isTRUE(observed_first) && !is.null(t_cond) && is.finite(t_cond)) {
        p <- p + ggplot2::annotate(
            "rect", xmin = -Inf, xmax = t_cond, ymin = -Inf, ymax = Inf,
            fill = "gray90", alpha = 0.3, colour = NA
        )
    }

    # Add credible bands (ribbons)
    if (ci_type %in% c("ribbon", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            alpha_level <- prediction_style$alpha / length(ci_levels)
            p_names <- .quantile_names_from_ci(level)

            if (all(p_names %in% names(df_pred))) {
                df_ribbon <- df_pred[, c("time", "marker", p_names), drop = FALSE]
                names(df_ribbon)[3:4] <- c("ymin", "ymax")

                p <- p + ggplot2::geom_ribbon(
                    ggplot2::aes(x = .data$time, ymin = .data$ymin, ymax = .data$ymax,
                                 fill = .data$marker),
                    data = df_ribbon,
                    alpha = alpha_level,
                    color = NA
                )
            }
        }
    }

    # Add quantile lines (if requested)
    if (ci_type %in% c("line", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            q_cols <- .quantile_names_from_ci(level)
            for (qcol in q_cols) {
                if (qcol %in% names(df_pred)) {
                    df_line <- df_pred[, c("time", "marker", qcol), drop = FALSE]
                    names(df_line)[3] <- "value"
                    p <- p + ggplot2::geom_line(
                        ggplot2::aes(x = .data$time, y = .data$value, color = .data$marker),
                        data = df_line,
                        linetype = "dashed",
                        alpha = 0.5,
                        linewidth = 0.4
                    )
                }
            }
        }
    }

    # Add median/mean trajectory
    median_col <- .quantile_name_from_prob(0.5)
    if (median_col %in% names(df_pred)) {
        p <- p + ggplot2::geom_line(
            ggplot2::aes(x = .data$time, y = .data[[median_col]], color = .data$marker),
            data = df_pred,
            linewidth = prediction_style$linewidth %||% 0.8
        )
    }

    # Add smoothed trajectory (optional)
    if (isTRUE(smooth_trajectory) && nrow(df_pred) > 3 && median_col %in% names(df_pred)) {
        if (smooth_method == "loess") {
            p <- p + ggplot2::geom_smooth(
                ggplot2::aes(x = .data$time, y = .data[[median_col]], color = .data$marker),
                data = df_pred,
                method = "loess",
                se = FALSE,
                span = smooth_span,
                alpha = 0.7,
                linetype = "dotted",
                linewidth = 0.6
            )
        }
    }

    # Add observed points and line
    if (isTRUE(show_data) && !is.null(df_obs) && nrow(df_obs) > 0) {
        if (isTRUE(show_observed_line)) {
            p <- p + ggplot2::geom_line(
                ggplot2::aes(x = .data$time, y = .data$value,
                             color = .data$marker, linetype = .data$type),
                data = df_obs,
                alpha = observed_style$alpha %||% 0.6,
                linewidth = 0.5
            ) +
                ggplot2::scale_linetype_manual(values = c(
                    "observed" = "solid",
                    "observed_future" = "dashed",
                    "fitted" = "solid",
                    "fitted_future" = "dashed"
                ))
        }

        p <- p + ggplot2::geom_point(
            ggplot2::aes(x = .data$time, y = .data$value, color = .data$marker),
            data = df_obs,
            shape = observed_style$shape %||% 21,
            size = observed_style$size %||% 2,
            alpha = observed_style$alpha %||% 0.6,
            fill = NA,
            stroke = 1
        )
    }

    # Faceting
    if (facet_by == "marker" && length(unique(df_pred$marker)) > 1) {
        p <- p + ggplot2::facet_wrap(~marker, scales = facet_scales)
    }

    # Color and fill scales
    if (!is.null(palette_marker)) {
        if (is.function(palette_marker)) {
            p <- p + palette_marker()
        } else if (is.character(palette_marker)) {
            unique_markers <- unique(c(df_pred$marker, if (!is.null(df_obs)) df_obs$marker else NULL))
            unique_markers <- unique_markers[!is.na(unique_markers)]
            if (length(unique_markers) > 0) {
                p <- p + ggplot2::scale_color_manual(values = palette_marker)
            }
        }
    } else {
        p <- p + ggplot2::scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)
    }

    scale_label <- switch(scale_long,
                          linpred = "linpred",
                          epred = "epred",
                          predict = "predict",
                          scale_long)

    y_label <- .latex_label_long(resp_var, scale_label)

    p <- p + ggplot2::scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
        ggplot2::labs(
            title = paste("Subject", id, "- Longitudinal Trajectory"),
            x = "Time",
            y = y_label,
            color = marker_var,
            fill = marker_var
        ) +
        theme_fn()

    p
}

.available_longitudinal_scales <- function(x) {
    scales <- x$metadata$scales %||% x$metadata$scale
    if (is.null(scales) && !is.null(x$quantiles$longitudinal) && "scale" %in% names(x$quantiles$longitudinal)) {
        scales <- unique(as.character(stats::na.omit(x$quantiles$longitudinal$scale)))
    }
    if (is.null(scales) && !is.null(x$predictions$longitudinal) && "scale" %in% names(x$predictions$longitudinal)) {
        scales <- unique(as.character(stats::na.omit(x$predictions$longitudinal$scale)))
    }
    scales <- unique(as.character(scales %||% "epred"))
    valid <- c("epred", "linpred", "predict")
    scales <- intersect(scales, valid)
    if (length(scales) == 0) {
        cli::cli_abort(c(
            x = "No valid longitudinal prediction scales are available in this object.",
            i = "Run {.fn predict} with {.arg process = 'longitudinal'} and a valid {.arg scale}."
        ))
    }
    scales
}

.available_plot_markers <- function(x) {
    marker_sources <- list(
        x$quantiles$longitudinal,
        x$predictions$longitudinal,
        x$quantiles$longitudinal_fitted,
        x$data$longitudinal
    )
    markers <- unique(unlist(lapply(marker_sources, function(df) {
        if (is.null(df) || !is.data.frame(df) || !"marker" %in% names(df)) return(character(0))
        as.character(stats::na.omit(df$marker))
    }), use.names = FALSE))
    markers[!is.na(markers) & nzchar(markers)]
}

.normalize_plot_marker_filter <- function(marker, available_markers = NULL) {
    if (is.null(marker) || length(marker) == 0 || all(is.na(marker))) {
        return(NULL)
    }
    marker <- unique(as.character(marker[!is.na(marker)]))
    if (length(marker) == 0) return(NULL)
    if (!is.null(available_markers) && length(available_markers) > 0) {
        missing_markers <- setdiff(marker, available_markers)
        if (length(missing_markers) > 0) {
            cli::cli_abort(c(
                x = "Unknown {.arg marker} level(s): {paste(missing_markers, collapse = ', ')}.",
                i = "Available markers: {paste(available_markers, collapse = ', ')}."
            ))
        }
    }
    marker
}

# ============================================================================
# Helper: Survival Plot for Single Subject
# ============================================================================
.plot_survival_single <- function(
    x, id,
    ci_levels, ci_type,
    prediction_style, theme_fn
) {

    # Extract quantile data
    quant_df <- x$quantiles$survival
    if (is.null(quant_df)) quant_df <- x$predictions$survival
    if (is.null(quant_df)) {
        cli::cli_warn(c(
            x = "No quantile or prediction data available for survival for subject {id}.",
            i = "Ensure survival predictions were requested."
        ))
        return(NULL)
    }

    quant_df <- quant_df[quant_df$id == id, , drop = FALSE]
    quant_df <- .ensure_plot_quantiles(quant_df, x$draws$survival[[as.character(id)]], id, ci_levels)
    if (nrow(quant_df) == 0) {
        cli::cli_warn(c(
            x = "No survival data found for subject {id}.",
            i = "Ensure survival predictions were requested."
        ))
        return(NULL)
    }

    # Initialise plot with data and base aesthetics (use q50 as y for survival)
    median_col <- .quantile_name_from_prob(0.5)
    if (median_col %in% names(quant_df)) {
        p <- ggplot2::ggplot(data = quant_df, ggplot2::aes(x = time, y = .data[[median_col]]))
    } else {
        p <- ggplot2::ggplot(data = quant_df, ggplot2::aes(x = time, y = mean))
    }

    # Add credible bands (ribbons)
    if (ci_type %in% c("ribbon", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            alpha_level <- prediction_style$alpha / length(ci_levels)
            p_names <- .quantile_names_from_ci(level)

            if (all(p_names %in% names(quant_df))) {
                p <- p + ggplot2::geom_ribbon(
                    ggplot2::aes(
                        x = .data$time,
                        ymin = .data[[p_names[1]]],
                        ymax = .data[[p_names[2]]]
                    ),
                    fill = prediction_style$fill %||% "steelblue",
                    alpha = alpha_level,
                    color = NA
                )
            }
        }
    }

    # Add median line (via base aes) with styling
    p <- p + ggplot2::geom_line(
        color = prediction_style$color %||% "steelblue",
        linewidth = prediction_style$linewidth %||% 0.8
    )

    # Add quantile lines
    if (ci_type %in% c("line", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            q_cols <- .quantile_names_from_ci(level)
            for (qcol in q_cols) {
                if (qcol %in% names(quant_df)) {
                    p <- p + ggplot2::geom_line(
                        ggplot2::aes(x = .data$time, y = .data[[qcol]]),
                        color = prediction_style$color %||% "steelblue",
                        linetype = "dashed",
                        alpha = 0.5,
                        linewidth = 0.4
                    )
                }
            }
        }
    }

    # Finalize
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        ggplot2::labs(
            title = paste("Subject", id, "- Conditional Survival"),
            x = "Time",
            y = .latex_label_surv()
        ) +
        theme_fn()

    p
}

# ============================================================================
# Helper: Cumulative Hazard Plot for Single Subject
# ============================================================================
.plot_cumhaz_single <- function(
    x, id,
    ci_levels, ci_type,
    prediction_style, theme_fn
) {

    quant_df <- x$quantiles$cumhaz
    if (is.null(quant_df)) quant_df <- x$predictions$cumhaz
    if (is.null(quant_df)) {
        cli::cli_warn(c(
            x = "No quantile or prediction data available for cumulative hazard for subject {id}.",
            i = "Ensure survival predictions were requested."
        ))
        return(NULL)
    }

    quant_df <- quant_df[quant_df$id == id, , drop = FALSE]
    quant_df <- .ensure_plot_quantiles(quant_df, x$draws$cumhaz[[as.character(id)]], id, ci_levels)
    if (nrow(quant_df) == 0) {
        cli::cli_warn(c(
            x = "No cumulative hazard data found for subject {id}.",
            i = "Ensure survival predictions were requested."
        ))
        return(NULL)
    }

    median_col <- .quantile_name_from_prob(0.5)
    if (median_col %in% names(quant_df)) {
        p <- ggplot2::ggplot(data = quant_df, ggplot2::aes(x = time, y = .data[[median_col]]))
    } else {
        p <- ggplot2::ggplot(data = quant_df, ggplot2::aes(x = time, y = mean))
    }

    if (ci_type %in% c("ribbon", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            alpha_level <- prediction_style$alpha / length(ci_levels)
            p_names <- .quantile_names_from_ci(level)
            if (all(p_names %in% names(quant_df))) {
                p <- p + ggplot2::geom_ribbon(
                    ggplot2::aes(
                        x = .data$time,
                        ymin = .data[[p_names[1]]],
                        ymax = .data[[p_names[2]]]
                    ),
                    fill = prediction_style$fill %||% "steelblue",
                    alpha = alpha_level,
                    color = NA
                )
            }
        }
    }

    p <- p + ggplot2::geom_line(
        color = prediction_style$color %||% "steelblue",
        linewidth = prediction_style$linewidth %||% 0.8
    )

    if (ci_type %in% c("line", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            q_cols <- .quantile_names_from_ci(level)
            for (qcol in q_cols) {
                if (qcol %in% names(quant_df)) {
                    p <- p + ggplot2::geom_line(
                        ggplot2::aes(x = .data$time, y = .data[[qcol]]),
                        color = prediction_style$color %||% "steelblue",
                        linetype = "dashed",
                        alpha = 0.5,
                        linewidth = 0.4
                    )
                }
            }
        }
    }

    p <- p + ggplot2::labs(
        title = paste("Subject", id, "- Conditional Cumulative Hazard"),
        x = "Time",
        y = "Cumulative hazard"
    ) + theme_fn()

    p
}

# ============================================================================
# Helper: resolve conditioning time per subject
# ============================================================================
.conditioning_time_for_id <- function(x, id) {
    tmap <- x$metadata$conditioning_time_by_id
    if (!is.null(tmap) && as.character(id) %in% names(tmap)) {
        return(as.numeric(tmap[[as.character(id)]]))
    }
    x$metadata$conditioning_time
}

# ============================================================================
# Helper: validate plot selection
# ============================================================================
.validate_plot_type <- function(type) {
    allowed <- c("cumhaz", "longitudinal", "survival")
    type <- unique(type)
    bad <- setdiff(type, allowed)
    if (length(bad) > 0) {
        cli::cli_abort(c(
            x = "Unknown plot type(s): {paste(bad, collapse = ', ')}.",
            i = "Use one or more of: {paste(allowed, collapse = ', ')}."
        ))
    }
    type
}

# ============================================================================
# Helper: ensure quantiles include requested CI levels
# ============================================================================
.ensure_plot_quantiles <- function(quant_df, draw_entry, id, ci_levels) {
    probs <- .quantile_probs_from_ci_plot(ci_levels)
    q_names <- .quantile_colnames(probs)
    missing_cols <- setdiff(q_names, names(quant_df))
    if (length(missing_cols) == 0) return(quant_df)

    if (is.null(draw_entry) || is.null(draw_entry$matrix)) {
        return(quant_df)
    }
    qdf <- .compute_quantiles_surv(draw_entry$matrix, draw_entry$time, id, probs = probs)
    if (is.null(quant_df) || nrow(quant_df) == 0) return(qdf)

    qdf <- qdf[, c("time", q_names), drop = FALSE]
    merged <- merge(quant_df, qdf, by = "time", all.x = TRUE, sort = FALSE)
    merged
}

.quantile_probs_from_ci_plot <- function(ci_levels) {
    bounds <- unlist(lapply(ci_levels, function(level) c((1 - level) / 2, (1 + level) / 2)))
    sort(unique(c(0.5, bounds)))
}

# ============================================================================
# Helper: combine plots using patchwork/cowplot
# ============================================================================
.combine_plot_grid <- function(plots, ncol = NULL, nrow = NULL, ..., fallback = c("flat", "input")) {
    fallback <- match.arg(fallback)
    plot_list <- .flatten_plot_list(plots)
    if (length(plot_list) == 0) return(NULL)

    if (requireNamespace("patchwork", quietly = TRUE)) {
        return(patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow, ...))
    }
    # if (requireNamespace("cowplot", quietly = TRUE)) {
    #     return(cowplot::plot_grid(plotlist = plot_list, ncol = ncol %||% length(plot_list)))
    # }
    # cli::cli_warn(c(
    #     x = "Neither {.pkg patchwork} nor {.pkg cowplot} is available for combined plots.",
    #     i = "Returning uncombined plots instead."
    # ))
    cli::cli_warn(c(
        x = "{.pkg patchwork} is not installed for combined plots.",
        i = "Returning uncombined plots instead."
    ))
    if (identical(fallback, "input")) plots else plot_list
}

.collect_plots_by_subject <- function(plots_list, which) {
    out <- list()
    for (id in names(plots_list)) {
        sub_plots <- plots_list[[id]]
        for (typ in which) {
            if (typ == "longitudinal") {
                long_keys <- names(sub_plots)[grepl("^longitudinal", names(sub_plots))]
                if (length(long_keys) == 0 && !is.null(sub_plots$longitudinal)) {
                    long_keys <- "longitudinal"
                }
                for (k in long_keys) {
                    out[[paste(id, k, sep = "_")]] <- sub_plots[[k]]
                }
            } else if (!is.null(sub_plots[[typ]])) {
                out[[paste(id, typ, sep = "_")]] <- sub_plots[[typ]]
            }
        }
    }
    out
}

.flatten_plot_list <- function(plots) {
    if (length(plots) == 0) return(list())
    if (all(vapply(plots, function(p) inherits(p, "gg"), logical(1)))) return(plots)
    out <- list()
    for (nm in names(plots)) {
        p <- plots[[nm]]
        if (inherits(p, "gg")) {
            out[[nm]] <- p
        } else if (is.list(p)) {
            inner <- .flatten_plot_list(p)
            for (k in names(inner)) out[[paste(nm, k, sep = "_")]] <- inner[[k]]
        }
    }
    out
}

# ============================================================================
# Helper: Combine longitudinal and survival plots
# ============================================================================
.combine_long_surv_plots <- function(p_long, p_surv) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        cli::cli_warn(c(
            x = "Package {.pkg patchwork} not available; returning list of plots instead.",
            i = "Install {.pkg patchwork} to enable combined plots."
        ))
        return(list(longitudinal = p_long, survival = p_surv))
    }

    p_long / p_surv
}

.validate_ci_levels_plot <- function(ci_levels, default_levels = NULL) {
    if (is.null(ci_levels) || length(ci_levels) == 0) {
        if (!is.null(default_levels)) {
            ci_levels <- default_levels
        } else {
            ci_levels <- c(0.5, 0.95)
        }
    }
    if (!is.numeric(ci_levels)) {
        cli::cli_abort(c(
            x = "{.arg ci_levels} must be numeric.",
            i = "Provide values between 0 and 1."
        ))
    }
    if (any(ci_levels <= 0 | ci_levels >= 1)) {
        cli::cli_abort(c(
            x = "{.arg ci_levels} must be strictly between 0 and 1.",
            i = "Example: ci_levels = c(0.5, 0.8, 0.95)."
        ))
    }
    sort(unique(ci_levels))
}

.quantile_name_from_prob <- function(prob) {
    pct <- 100 * prob
    lbl <- formatC(pct, format = "fg", digits = 6)
    lbl <- gsub("\\s+", "", lbl)
    lbl <- sub("\\.?0+$", "", lbl)
    paste0("q", lbl)
}

.quantile_names_from_ci <- function(level) {
    .quantile_name_from_prob(c((1 - level) / 2, (1 + level) / 2))
}

.latex_label_long <- function(resp_var, scale_label) {
    base_label <- paste0(resp_var, " (", scale_label, ")")
    if (!is.na(resp_var) && requireNamespace("latex2exp", quietly = TRUE)) {
        latex2exp::TeX(paste0("$", resp_var, "$\\,(", scale_label, ")"))
    } else {
        base_label
    }
}

.latex_label_surv <- function() {
    if (requireNamespace("latex2exp", quietly = TRUE)) {
        latex2exp::TeX("$S(t \\mid T_{cond})$")
    } else {
        "S(t | T_cond)"
    }
}

# ============================================================================
# Diagnostic plots for JoinMeFit
# ============================================================================

#' Plot diagnostics, fitted trajectories, and association curves for JoinMe fits
#'
#' @param x A fitted object of class `JoinMeFit`.
#' @param type Plot type. Diagnostic types are `"rhat"`, `"ess_bulk"`,
#'   `"ess_tail"`, `"mcse_mean"`, `"mcse_sd"`, `"running_mean"`, and
#'   `"running_quantile"`. Fitted-data types are `"longitudinal"`,
#'   `"survival"`, `"cumhaz"`, and `"association"`. The compatibility alias
#'   `"longitudinal_heatmap"` is treated as
#'   `type = "longitudinal", longitudinal_style = "heatmap"`.
#' @param pars Optional character vector of parameter names to include for
#'   diagnostic plots.
#' @param regex_pars Optional regular expression for parameter selection for
#'   diagnostic plots.
#' @param draws Optional number of posterior draws to subset for diagnostics.
#' @param seed Random seed for draw subsetting and fitted plotting.
#' @param max_vars Maximum number of parameters for running diagnostics.
#' @param quantile_probs Numeric vector of quantile probabilities for running
#'   quantile plots.
#' @param subject Optional vector of subject ids for fitted-data plots.
#' @param marker Optional marker subset for longitudinal and association plots.
#'   Use a specific marker level to plot only that marker, or `NA` to plot all
#'   markers.
#' @param scale Longitudinal scale for fitted-data plots. One of `"epred"`,
#'   `"linpred"`, or `"predict"`.
#' @param longitudinal_style Display style for fitted longitudinal plots.
#'   `"curves"` shows the existing separated trajectory curves with credible
#'   bands. `"heatmap"` shows marker-level mean change over time using draw-level
#'   aggregation across the selected subjects. The heatmap is therefore an
#'   alternative display for the longitudinal process rather than a separate
#'   modelling target.
#' @param condition Optional conditioning specification for fitted-data plots.
#'   Supply either a named list for one condition or a data frame with one row
#'   per condition, such as the output of [make_conditions()]. Each condition
#'   row may contain baseline covariates from the longitudinal data, the event
#'   data, or both. For ordinary fitted trajectories, multiple rows return a
#'   named list of condition-specific plots unless a single plot type is
#'   requested with `combined = TRUE`, in which case the condition-specific
#'   panels are combined. For `longitudinal_style = "heatmap"`, multiple
#'   condition rows are shown as separate facets.
#' @param conditioning Conditioning-time convention when `time_start` is not
#'   supplied. `"last"` uses each subject's last observed longitudinal time,
#'   `"origin"` uses the first observed time, and `"auto"` chooses `"last"`
#'   when longitudinal curves are requested and `"origin"` otherwise.
#' @param time_start Optional conditioning time passed to fitted-data
#'   predictions. Can be scalar numeric, named numeric vector, or an event-data
#'   column name.
#' @param times Optional time grid for fitted-data prediction.
#' @param time_horizon Optional prediction horizon for fitted-data prediction.
#' @param pred_control Named list passed to `predict.JoinMeFit()` for fitted-data
#'   plotting. By default this uses a modest draw count for plotting.
#' @param threshold Posterior sign-certainty threshold for
#'   `longitudinal_style = "heatmap"`. A tile is shown as significant when the
#'   posterior probability of either a positive or a negative change, relative
#'   to the earliest plotted time for that marker, is at least `1 - threshold`.
#'   Tiles that do not meet that criterion are drawn transparently.
#' @param smooth_trajectory,smooth_method,smooth_span,ci_levels,ci_type,
#'   observed_first,facet_by,facet_scales,combined,show_data,
#'   show_observed_line,observed_style,prediction_style,theme_fn,
#'   palette_marker Passed to `plot.JoinMeDynPred()` for fitted-data plots.
#' @param association_options Named list of options for
#'   `type = "association"`. Supported entries are
#'   `association_term`, `association_grid`, `association_range`,
#'   `association_points`, and `association_metric`. When
#'   `association_grid` is not supplied, `plot.JoinMeFit()` uses cached
#'   model-implied raw support from the fitted association channels and only
#'   falls back to knot support or observed-data heuristics when that cache is
#'   unavailable. `association_range`, when supplied, overrides those default
#'   range heuristics and constructs an evenly spaced raw-scale grid over the
#'   requested interval.
#'   `association_metric = "hazard"` plots the posterior contribution used by
#'   the fitted model. For covariance-style channels (`corr`, `vcov`) this
#'   contribution is zero-referenced at raw value `0` before multiplying by
#'   \eqn{\alpha}; `association_metric = "transform"` plots the transform $f(x)$
#'   alone and therefore omits the association-coefficient sign.
#' @param ... Unused.
#'
#' @return A `ggplot` object, a combined plot, or a named list of plots.
#' @export
plot.JoinMeFit <- function(x,
                           type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean", "running_quantile"),
                           pars = NULL, regex_pars = NULL, draws = NULL, seed = 1, max_vars = 4,
                           quantile_probs = c(0.1, 0.5, 0.9),
                           subject = NULL,
                           marker = NA,
                           scale = NULL,
                           longitudinal_style = c("curves", "heatmap"),
                           condition = NULL,
                           conditioning = c("auto", "last", "origin"),
                           time_start = NULL,
                           times = NULL,
                           time_horizon = NULL,
                           pred_control = list(n_samples = 100),
                           threshold = 0.05,
                           smooth_trajectory = TRUE,
                           smooth_method = c("loess", "spline"),
                           smooth_span = 0.3,
                           ci_levels = c(0.5, 0.95),
                           ci_type = c("ribbon", "line", "both"),
                           observed_first = TRUE,
                           facet_by = c("marker", "none"),
                           facet_scales = "free_y",
                           combined = TRUE,
                           show_data = TRUE,
                           show_observed_line = TRUE,
                           observed_style = list(
                               color = "black",
                               shape = 21,
                               size = 2,
                               alpha = 0.6
                           ),
                           prediction_style = list(
                               color = "steelblue",
                               fill = "steelblue",
                               linewidth = 0.8,
                               alpha = 0.2
                           ),
                           theme_fn = ggplot2::theme_bw,
                           palette_marker = NULL,
                           association_options = list(),
                           ...) {
    dots <- list(...)
    diagnostic_types <- c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean", "running_quantile")
    fitted_types <- c("longitudinal", "survival", "cumhaz", "association", "longitudinal_heatmap")

    if (!inherits(x, "JoinMeFit")) {
        cli::cli_abort("{.arg x} must be a JoinMeFit object.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }

    if ("which" %in% names(dots)) {
        type <- dots$which
    }

    type <- tryCatch(match.arg(unique(as.character(type)), c(diagnostic_types, fitted_types), several.ok = TRUE), error = function(e) {
        cli::cli_abort(c(
            x = "Unknown {.arg type}: {paste(type, collapse = ', ')}.",
            i = "Use one or more of: {paste(c(diagnostic_types, fitted_types), collapse = ', ')}."
        ))
    })

    if (any(type %in% diagnostic_types) && any(type %in% fitted_types)) {
        cli::cli_abort(c(
            x = "Diagnostic and fitted-data plot types cannot be mixed in one call.",
            i = "Call {.fn plot} separately for diagnostics and fitted trajectories."
        ))
    }
    if (all(type %in% diagnostic_types)) {
        if (length(type) != 1L) {
            cli::cli_abort(c(
                x = "Diagnostic plotting expects a single {.arg type}.",
                i = "Example: plot(fit, type = 'rhat')."
            ))
        }
        return(.plot_joinmefit_diagnostic(
            x = x,
            type = type,
            pars = pars,
            regex_pars = regex_pars,
            draws = draws,
            seed = seed,
            max_vars = max_vars,
            quantile_probs = quantile_probs
        ))
    }

    longitudinal_style <- match.arg(longitudinal_style)
    conditioning <- match.arg(conditioning)
    smooth_method <- match.arg(smooth_method)
    facet_by <- match.arg(facet_by)
    ci_type <- match.arg(ci_type)

    if ("longitudinal_heatmap" %in% type) {
        longitudinal_style <- "heatmap"
        type <- unique(c(setdiff(type, "longitudinal_heatmap"), "longitudinal"))
    }

    if (identical(longitudinal_style, "heatmap") && (!all(type %in% c("longitudinal")) || length(type) != 1L)) {
        cli::cli_abort(c(
            x = "{.arg longitudinal_style = 'heatmap'} requires {.arg type = 'longitudinal'}.",
            i = "Use the heatmap as an alternative longitudinal display, not in combination with survival, cumulative hazard, or association plots."
        ))
    }

    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    marker_levels <- if (!is.null(x$dataLong) && marker_var %in% names(x$dataLong)) {
        unique(as.character(stats::na.omit(x$dataLong[[marker_var]])))
    } else {
        character(0)
    }
    marker <- .normalize_plot_marker_filter(marker, available_markers = marker_levels)

    if (identical(type, "association") || (length(type) == 1L && type[[1]] == "association")) {
        association_opts <- .normalize_joinmefit_association_options(
            association_options = association_options,
            dots = dots
        )
        return(.plot_joinmefit_association(
            x = x,
            association_term = association_opts$association_term,
            marker = marker,
            association_grid = association_opts$association_grid,
            association_range = association_opts$association_range,
            association_points = association_opts$association_points,
            ci_levels = ci_levels,
            ci_type = ci_type,
            association_metric = association_opts$association_metric,
            prediction_style = prediction_style,
            theme_fn = theme_fn,
            combined = combined,
            seed = seed
        ))
    }

    .plot_joinmefit_fitted(
        x = x,
        type = type,
        subject = subject,
        marker = marker,
        scale = scale,
        longitudinal_style = longitudinal_style,
        condition = condition,
        conditioning = conditioning,
        time_start = time_start,
        times = times,
        time_horizon = time_horizon,
        pred_control = pred_control,
        threshold = threshold,
        seed = seed,
        smooth_trajectory = smooth_trajectory,
        smooth_method = smooth_method,
        smooth_span = smooth_span,
        ci_levels = ci_levels,
        ci_type = ci_type,
        observed_first = observed_first,
        facet_by = facet_by,
        facet_scales = facet_scales,
        combined = combined,
        show_data = show_data,
        show_observed_line = show_observed_line,
        observed_style = observed_style,
        prediction_style = prediction_style,
        theme_fn = theme_fn,
        palette_marker = palette_marker
    )
}

.normalize_joinmefit_plot_conditions <- function(condition) {
    if (is.null(condition)) {
        return(list(rows = list(NULL), labels = NULL))
    }

    if (inherits(condition, "data.frame")) {
        condition_df <- condition
    } else if (is.list(condition)) {
        condition_df <- as.data.frame(condition, stringsAsFactors = FALSE)
    } else {
        cli::cli_abort(c(
            x = "{.arg condition} must be NULL, a named list, or a data frame.",
            i = "Use a one-row named list for one condition, or pass the result of {.fn make_conditions}."
        ))
    }

    if (!nrow(condition_df)) {
        cli::cli_abort(c(
            x = "{.arg condition} must contain at least one row.",
            i = "Provide one or more condition rows to define the plotted covariate profile."
        ))
    }

    labels <- rownames(condition_df)
    labels[is.na(labels) | !nzchar(labels)] <- paste("Condition", which(is.na(labels) | !nzchar(labels)))
    if (is.null(labels)) {
        labels <- paste("Condition", seq_len(nrow(condition_df)))
    }

    rows <- lapply(seq_len(nrow(condition_df)), function(i) {
        row_i <- condition_df[i, , drop = FALSE]
        keep <- vapply(row_i, function(col) !all(is.na(col)), logical(1))
        as.list(row_i[, keep, drop = FALSE])
    })

    list(rows = rows, labels = labels)
}

.coerce_joinmefit_condition_value <- function(template, value, n) {
    value <- value[[1L]]

    if (is.factor(template)) {
        return(factor(rep(as.character(value), n), levels = levels(template), ordered = is.ordered(template)))
    }
    if (inherits(template, "Date")) {
        return(rep(as.Date(value), n))
    }
    if (inherits(template, "POSIXct")) {
        tz <- attr(template, "tzone") %||% ""
        return(rep(as.POSIXct(value, tz = tz), n))
    }
    if (is.integer(template)) {
        return(as.integer(rep(value, n)))
    }
    if (is.numeric(template)) {
        return(as.numeric(rep(value, n)))
    }
    if (is.logical(template)) {
        return(as.logical(rep(value, n)))
    }

    rep(value, n)
}

.apply_joinmefit_plot_condition <- function(data_long, data_event, condition_row, protected_columns = character(0)) {
    if (is.null(condition_row) || !length(condition_row)) {
        return(list(longitudinal = data_long, event = data_event))
    }

    condition_names <- names(condition_row)
    allowed_long <- setdiff(names(data_long), protected_columns)
    allowed_event <- setdiff(names(data_event), protected_columns)
    unknown <- setdiff(condition_names, union(allowed_long, allowed_event))
    if (length(unknown) > 0L) {
        cli::cli_abort(c(
            x = "Unknown {.arg condition} column{?s}: {paste(unknown, collapse = ', ')}.",
            i = "Condition columns must match baseline covariates in the longitudinal or event data."
        ))
    }

    for (nm in condition_names) {
        if (nm %in% allowed_long && nrow(data_long) > 0L) {
            data_long[[nm]] <- .coerce_joinmefit_condition_value(data_long[[nm]], condition_row[[nm]], nrow(data_long))
        }
        if (nm %in% allowed_event && nrow(data_event) > 0L) {
            data_event[[nm]] <- .coerce_joinmefit_condition_value(data_event[[nm]], condition_row[[nm]], nrow(data_event))
        }
    }

    list(longitudinal = data_long, event = data_event)
}

.joinmefit_longitudinal_draw_scale <- function(draw_obj, scale) {
    if (is.null(draw_obj)) {
        return(NULL)
    }
    if (!is.null(draw_obj$matrix) && !is.null(draw_obj$scale)) {
        return(draw_obj)
    }
    draw_obj[[scale]] %||% NULL
}

.joinmefit_longitudinal_heatmap_data <- function(pred, scale, threshold, marker = NULL) {
    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1) {
        cli::cli_abort(c(
            x = "{.arg threshold} must be a single number strictly between 0 and 1.",
            i = "A common choice is {.code threshold = 0.05} for a 95% posterior sign-certainty rule."
        ))
    }

    draws_long <- pred$draws$longitudinal
    if (is.null(draws_long) || !length(draws_long)) {
        cli::cli_abort(c(
            x = "Longitudinal draw-level predictions are required for {.val longitudinal_heatmap}.",
            i = "Call {.fn plot} on a fitted JoinMe model so the helper can build draw-level predictions automatically."
        ))
    }

    marker_var <- pred$metadata$marker_var %||% "marker"
    marker_source <- pred$data$longitudinal[[marker_var]]
    marker_levels <- if (is.factor(marker_source)) {
        levels(marker_source)
    } else {
        sort(unique(as.character(stats::na.omit(marker_source))))
    }

    reference_layout <- NULL
    aggregated_matrix <- NULL
    n_subjects <- 0L

    for (id in names(draws_long)) {
        draw_scale <- .joinmefit_longitudinal_draw_scale(draws_long[[id]], scale)
        if (is.null(draw_scale) || is.null(draw_scale$matrix)) {
            next
        }

        layout_i <- data.frame(
            time = as.numeric(draw_scale$time),
            marker = marker_levels[as.integer(draw_scale$marker_idx)],
            stringsAsFactors = FALSE
        )

        if (!is.null(marker)) {
            keep_cols <- as.character(layout_i$marker) %in% marker
            layout_i <- layout_i[keep_cols, , drop = FALSE]
            draw_scale$matrix <- draw_scale$matrix[, keep_cols, drop = FALSE]
        }

        if (!nrow(layout_i)) {
            next
        }

        if (is.null(reference_layout)) {
            reference_layout <- layout_i
            aggregated_matrix <- draw_scale$matrix
        } else {
            same_marker <- identical(as.character(reference_layout$marker), as.character(layout_i$marker))
            same_time <- isTRUE(all.equal(reference_layout$time, layout_i$time, tolerance = 1e-10))
            if (!same_marker || !same_time) {
                cli::cli_abort(c(
                    x = "All subjects must share the same time-marker prediction grid for {.val longitudinal_heatmap}.",
                    i = "Supply a common {.arg times} grid or use subjects with compatible marker schedules."
                ))
            }
            aggregated_matrix <- aggregated_matrix + draw_scale$matrix
        }
        n_subjects <- n_subjects + 1L
    }

    if (is.null(reference_layout) || is.null(aggregated_matrix) || n_subjects == 0L) {
        cli::cli_abort(c(
            x = "No longitudinal predictions were available for the requested heatmap.",
            i = "Check the requested {.arg subject}, {.arg marker}, and {.arg scale} filters."
        ))
    }

    aggregated_matrix <- aggregated_matrix / n_subjects
    by_marker <- split(seq_len(nrow(reference_layout)), as.character(reference_layout$marker))

    heatmap_df <- do.call(rbind, lapply(names(by_marker), function(marker_name) {
        idx <- by_marker[[marker_name]]
        idx <- idx[order(reference_layout$time[idx])]
        baseline <- aggregated_matrix[, idx[1L], drop = FALSE]
        change_matrix <- sweep(aggregated_matrix[, idx, drop = FALSE], 1L, baseline[, 1L], FUN = "-")
        prob_positive <- colMeans(change_matrix > 0, na.rm = TRUE)
        prob_negative <- colMeans(change_matrix < 0, na.rm = TRUE)
        posterior_sign_prob <- pmax(prob_positive, prob_negative)

        data.frame(
            time = reference_layout$time[idx],
            marker = marker_name,
            change = colMeans(change_matrix, na.rm = TRUE),
            posterior_sign_prob = posterior_sign_prob,
            significant = posterior_sign_prob >= (1 - threshold),
            alpha = ifelse(posterior_sign_prob >= (1 - threshold), 1, 0),
            stringsAsFactors = FALSE
        )
    }))

    marker_order <- data.frame(
        marker = names(split(heatmap_df$change, heatmap_df$marker)),
        order_value = vapply(split(abs(heatmap_df$change), heatmap_df$marker), max, numeric(1), na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    marker_order <- marker_order[order(-marker_order$order_value, marker_order$marker), , drop = FALSE]

    heatmap_df$marker <- factor(heatmap_df$marker, levels = marker_order$marker)
    heatmap_df
}

.plot_joinmefit_longitudinal_heatmap <- function(x, subject, marker, scale, condition, conditioning,
                                                 time_start, times, time_horizon, pred_control,
                                                 threshold, seed, ci_levels, theme_fn) {
    scale_use <- .normalize_prediction_scales(scale %||% "epred")[[1L]]
    condition_spec <- .normalize_joinmefit_plot_conditions(condition)

    plot_data <- do.call(rbind, lapply(seq_along(condition_spec$rows), function(i) {
        pred <- .build_joinmefit_plot_prediction(
            x = x,
            which = "longitudinal",
            subject = subject,
            scale = scale_use,
            condition = condition_spec$rows[[i]],
            conditioning = conditioning,
            time_start = time_start,
            times = times,
            time_horizon = time_horizon,
            pred_control = pred_control,
            seed = seed,
            ci_levels = ci_levels,
            pred_type = "per_marker_id"
        )

        data_i <- .joinmefit_longitudinal_heatmap_data(
            pred = pred,
            scale = scale_use,
            threshold = threshold,
            marker = marker
        )
        data_i$condition_label <- condition_spec$labels[[i]] %||% "Observed data"
        data_i
    }))

    p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = .data$time, y = .data$marker, fill = .data$change, alpha = .data$alpha)
    ) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::scale_alpha_identity() +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC",
            mid = "#F7F7F7",
            high = "#B2182B",
            midpoint = 0,
            name = paste0("Change (", scale_use, ")")
        ) +
        ggplot2::labs(
            x = "Time",
            y = "Marker",
            title = "Longitudinal mean change",
            subtitle = "Change is measured relative to the earliest plotted time for each marker"
        ) +
        theme_fn()

    if (length(unique(plot_data$condition_label)) > 1L) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$condition_label), scales = "free_y")
    }

    p
}

.normalize_joinmefit_association_options <- function(association_options = NULL, dots = list()) {
    if (is.null(association_options)) {
        association_options <- list()
    }
    if (!is.list(association_options)) {
        cli::cli_abort("{.arg association_options} must be a named list.")
    }

    if (length(association_options) > 0L) {
        nm <- names(association_options)
        if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
            cli::cli_abort("{.arg association_options} must be a named list.")
        }
    }

    alias_map <- c(
        term = "association_term",
        grid = "association_grid",
        range = "association_range",
        points = "association_points",
        metric = "association_metric"
    )
    valid_names <- c(unname(alias_map), names(alias_map))
    unknown <- setdiff(names(association_options), valid_names)
    if (length(unknown) > 0L) {
        cli::cli_abort(c(
            x = "Unknown {.arg association_options} entr{?y/ies}: {paste(unknown, collapse = ', ')}.",
            i = "Supported names are association_term, association_grid, association_range, association_points, association_metric, and their aliases term, grid, range, points, metric."
        ))
    }

    canonical_names <- ifelse(names(association_options) %in% names(alias_map), alias_map[names(association_options)], names(association_options))
    duplicated_names <- unique(canonical_names[duplicated(canonical_names)])
    if (length(duplicated_names) > 0L) {
        cli::cli_abort("{.arg association_options} cannot supply the same option more than once.")
    }
    names(association_options) <- canonical_names

    opts <- list(
        association_term = NULL,
        association_grid = NULL,
        association_range = NULL,
        association_points = 200,
        association_metric = "hazard"
    )
    opts[names(association_options)] <- association_options

    opts$association_metric <- match.arg(as.character(opts$association_metric)[1], c("hazard", "transform"))

    if (!is.null(opts$association_grid)) {
        opts$association_grid <- as.numeric(opts$association_grid)
        if (!length(opts$association_grid) || any(!is.finite(opts$association_grid))) {
            cli::cli_abort("{.arg association_grid} must be a numeric vector with finite values.")
        }
    }

    if (!is.null(opts$association_range)) {
        opts$association_range <- as.numeric(opts$association_range)
        if (length(opts$association_range) != 2L || any(!is.finite(opts$association_range))) {
            cli::cli_abort("{.arg association_range} must be a finite numeric vector of length 2.")
        }
        opts$association_range <- sort(opts$association_range)
        if (diff(opts$association_range) <= 0) {
            cli::cli_abort("{.arg association_range} must span two distinct values.")
        }
    }

    if (!is.null(opts$association_grid) && !is.null(opts$association_range)) {
        cli::cli_abort(c(
            x = "{.arg association_grid} and {.arg association_range} cannot be supplied together.",
            i = "Use {.arg association_grid} for an explicit evaluation grid, or {.arg association_range} with {.arg association_points} for an evenly spaced grid."
        ))
    }

    opts$association_points <- as.integer(opts$association_points)[1]
    if (!is.finite(opts$association_points) || opts$association_points < 2L) {
        cli::cli_abort("{.arg association_points} must be an integer greater than or equal to 2.")
    }

    opts
}

.plot_joinmefit_diagnostic <- function(x, type, pars, regex_pars, draws, seed, max_vars, quantile_probs) {
    if (type %in% c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd")) {
        df <- switch(type,
            rhat = stan_rhat.JoinMeFit(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed),
            ess_bulk = stan_ess.JoinMeFit(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "bulk"),
            ess_tail = stan_ess.JoinMeFit(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "tail"),
            mcse_mean = stan_mcse.JoinMeFit(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "mean"),
            mcse_sd = stan_mcse.JoinMeFit(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "sd")
        )
        if (nrow(df) == 0) {
            cli::cli_abort("No parameters found for diagnostic plotting.")
        }
        metric_name <- setdiff(names(df), "variable")
        names(df)[names(df) == metric_name] <- "metric"
        p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(.data$variable, .data$metric), y = .data$metric)) +
            ggplot2::geom_point(size = 1.6, alpha = 0.7, color = "steelblue") +
            ggplot2::coord_flip() +
            ggplot2::labs(x = NULL, y = type, title = paste("JoinMe diagnostics:", type)) +
            ggplot2::theme_minimal()
        if (type == "rhat") {
            p <- p + ggplot2::geom_hline(yintercept = 1.01, linetype = "dashed", color = "firebrick")
        }
        return(p)
    }

    vars <- posterior::variables(.get_draws_obj(x$fit))
    vars <- .filter_diag_vars(vars, pars, regex_pars)
    if (length(vars) == 0) {
        cli::cli_abort("No parameters found for running diagnostics.")
    }
    if (length(vars) > max_vars) {
        vars <- vars[seq_len(max_vars)]
        cli::cli_warn("Limiting running diagnostics to the first {max_vars} parameters.")
    }

    draws_obj <- .get_draws_obj(x$fit, variables = vars, draws = draws, seed = seed)
    arr <- posterior::as_draws_array(draws_obj)
    if (type == "running_mean") {
        df <- .running_mean_df(arr)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = factor(.data$chain))) +
            ggplot2::geom_line(alpha = 0.7) +
            ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
            ggplot2::labs(x = "Iteration", y = "Running mean", color = "Chain", title = "Running mean by chain") +
            ggplot2::theme_minimal()
        return(p)
    }

    df <- .running_quantile_df(arr, probs = quantile_probs)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = factor(.data$chain))) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::facet_grid(rows = ggplot2::vars(.data$variable), cols = ggplot2::vars(.data$stat), scales = "free_y") +
        ggplot2::labs(x = "Iteration", y = "Running quantile", color = "Chain", title = "Running quantiles by chain") +
        ggplot2::theme_minimal()
}

.plot_joinmefit_fitted <- function(x, type, subject, marker, scale, longitudinal_style, condition, conditioning, time_start, times, time_horizon,
                                   pred_control, threshold, seed, smooth_trajectory, smooth_method, smooth_span,
                                   ci_levels, ci_type, observed_first, facet_by, facet_scales, combined,
                                   show_data, show_observed_line, observed_style, prediction_style,
                                   theme_fn, palette_marker) {
    if (identical(longitudinal_style, "heatmap")) {
        return(.plot_joinmefit_longitudinal_heatmap(
            x = x,
            subject = subject,
            marker = marker,
            scale = scale,
            condition = condition,
            conditioning = conditioning,
            time_start = time_start,
            times = times,
            time_horizon = time_horizon,
            pred_control = pred_control,
            threshold = threshold,
            seed = seed,
            ci_levels = ci_levels,
            theme_fn = theme_fn
        ))
    }

    condition_spec <- .normalize_joinmefit_plot_conditions(condition)
    plot_types <- intersect(type, c("longitudinal", "survival", "cumhaz"))

    if (length(condition_spec$rows) > 1L) {
        plots <- lapply(seq_along(condition_spec$rows), function(i) {
            pred <- .build_joinmefit_plot_prediction(
                x = x,
                which = type,
                subject = subject,
                scale = scale,
                condition = condition_spec$rows[[i]],
                conditioning = conditioning,
                time_start = time_start,
                times = times,
                time_horizon = time_horizon,
                pred_control = pred_control,
                seed = seed,
                ci_levels = ci_levels,
                pred_type = "per_marker_id"
            )

            out_i <- plot(pred,
                type = plot_types,
                subject = subject,
                marker = if (is.null(marker)) NA else marker,
                scale = scale,
                smooth_trajectory = smooth_trajectory,
                smooth_method = smooth_method,
                smooth_span = smooth_span,
                ci_levels = ci_levels,
                ci_type = ci_type,
                observed_first = observed_first,
                facet_by = facet_by,
                facet_scales = facet_scales,
                combined = combined,
                show_data = show_data,
                show_observed_line = show_observed_line,
                observed_style = observed_style,
                prediction_style = prediction_style,
                theme_fn = theme_fn,
                palette_marker = palette_marker)

            if (inherits(out_i, "ggplot")) {
                out_i <- out_i + ggplot2::labs(subtitle = condition_spec$labels[[i]])
            }
            out_i
        })
        names(plots) <- condition_spec$labels

        if (isTRUE(combined) && length(plot_types) == 1L && all(vapply(plots, inherits, logical(1), what = "ggplot"))) {
            return(.combine_plot_grid(plots, fallback = "input"))
        }
        return(plots)
    }

    pred <- .build_joinmefit_plot_prediction(
        x = x,
        which = type,
        subject = subject,
        scale = scale,
        condition = condition_spec$rows[[1L]],
        conditioning = conditioning,
        time_start = time_start,
        times = times,
        time_horizon = time_horizon,
        pred_control = pred_control,
        seed = seed,
        ci_levels = ci_levels,
        pred_type = "per_marker_id"
    )

    plot(pred,
         type = plot_types,
         subject = subject,
         marker = if (is.null(marker)) NA else marker,
         scale = scale,
         smooth_trajectory = smooth_trajectory,
         smooth_method = smooth_method,
         smooth_span = smooth_span,
         ci_levels = ci_levels,
         ci_type = ci_type,
         observed_first = observed_first,
         facet_by = facet_by,
         facet_scales = facet_scales,
         combined = combined,
         show_data = show_data,
         show_observed_line = show_observed_line,
         observed_style = observed_style,
         prediction_style = prediction_style,
         theme_fn = theme_fn,
         palette_marker = palette_marker)
}

.build_joinmefit_plot_prediction <- function(x, which, subject, scale, condition, conditioning, time_start, times,
                                             time_horizon, pred_control, seed, ci_levels, pred_type = "per_marker_id") {
    id_var <- .joinmefit_call_arg_chr(x$call, "id_var", "id")
    time_var <- .joinmefit_call_arg_chr(x$call, "time_var", "time")
    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    response_var <- tryCatch(all.vars(x$formulaLong)[1], error = function(e) NA_character_)
    event_outcome_vars <- tryCatch(all.vars(x$formulaEvent[[2]]), error = function(e) character(0))

    data_long <- x$dataLong
    data_event <- x$dataEvent
    if (!is.null(subject)) {
        keep_ids <- as.character(subject)
        data_long <- data_long[as.character(data_long[[id_var]]) %in% keep_ids, , drop = FALSE]
        data_event <- data_event[as.character(data_event[[id_var]]) %in% keep_ids, , drop = FALSE]
    }

    data_long_plot <- data_long
    data_event_plot <- data_event

    protected_columns <- unique(stats::na.omit(c(id_var, time_var, marker_var, response_var, event_outcome_vars)))
    conditioned_data <- .apply_joinmefit_plot_condition(
        data_long = data_long,
        data_event = data_event,
        condition_row = condition,
        protected_columns = protected_columns
    )
    data_long <- conditioned_data$longitudinal
    data_event <- conditioned_data$event

    time_start_use <- .resolve_joinmefit_plot_time_start(
        data_long = data_long,
        id_var = id_var,
        time_var = time_var,
        which = which,
        conditioning = conditioning,
        time_start = time_start
    )

    process <- c(
        if (any(c("longitudinal", "longitudinal_heatmap") %in% which)) "longitudinal",
        if (any(c("longitudinal", "longitudinal_heatmap") %in% which) || any(c("survival", "cumhaz") %in% which)) "event"
    )
    process <- unique(process)

    scale_use <- if (any(c("longitudinal", "longitudinal_heatmap") %in% which)) .normalize_prediction_scales(scale %||% c("epred", "linpred", "predict")) else NULL
    control_use <- utils::modifyList(list(n_samples = 100L), pred_control %||% list())

    pred <- predict.JoinMeFit(
        object = x,
        newdataLong = data_long,
        newdataEvent = data_event,
        process = process,
        pred_type = pred_type,
        scale = scale_use,
        times = times,
        time_start = time_start_use,
        time_horizon = time_horizon,
        ci_levels = ci_levels,
        control = control_use,
        seed = seed
    )

    pred$data$longitudinal <- data_long_plot
    pred$data$event <- data_event_plot
    pred$metadata$plot_condition <- condition
    pred
}

.joinmefit_call_arg_chr <- function(call_obj, arg, default = NULL) {
    expr <- call_obj[[arg]]
    if (is.null(expr)) return(default)
    if (is.character(expr)) return(expr[[1]])
    txt <- trimws(paste(deparse(expr), collapse = ""))
    if (identical(txt, "") || identical(txt, "NULL")) default else txt
}

.resolve_joinmefit_plot_time_start <- function(data_long, id_var, time_var, which, conditioning, time_start) {
    if (!is.null(time_start)) return(time_start)
    if (!nrow(data_long)) return(NULL)

    conditioning <- match.arg(conditioning, c("auto", "last", "origin"))
    if (identical(conditioning, "auto")) {
        conditioning <- if ("longitudinal" %in% which) "last" else "origin"
    }

    split_time <- split(as.numeric(data_long[[time_var]]), as.character(data_long[[id_var]]))
    reducer <- if (identical(conditioning, "last")) max else min
    stats::setNames(
        vapply(split_time, function(tt) reducer(tt, na.rm = TRUE), numeric(1)),
        names(split_time)
    )
}

.plot_joinmefit_association <- function(x, association_term = NULL, marker = NULL, association_grid = NULL, association_range = NULL, association_points = 200,
                                        ci_levels = c(0.5, 0.95), ci_type = c("ribbon", "line", "both"),
                                        association_metric = c("hazard", "transform"),
                                        prediction_style = list(color = "steelblue", fill = "steelblue", linewidth = 0.8, alpha = 0.2),
                                        theme_fn = ggplot2::theme_bw, combined = TRUE, seed = 1) {
    ci_type <- match.arg(ci_type)
    association_metric <- match.arg(association_metric)
    assoc_terms <- .available_association_terms(x)
    if (!is.null(association_term)) {
        assoc_terms <- .expand_association_terms(association_term, assoc_terms)
    }
    if (length(assoc_terms) == 0) {
        cli::cli_abort(c(
            x = "No association terms available for plotting.",
            i = "Fit a model with association terms such as {.val cv_total} or {.val corr}."
        ))
    }

    plots <- lapply(assoc_terms, function(term) {
        .plot_joinmefit_association_single(
            x = x,
            term = term,
            marker = marker,
            association_grid = association_grid,
            association_range = association_range,
            association_points = association_points,
            ci_levels = ci_levels,
            ci_type = ci_type,
            association_metric = association_metric,
            prediction_style = prediction_style,
            theme_fn = theme_fn,
            seed = seed
        )
    })
    names(plots) <- assoc_terms

    if (length(plots) == 1L) return(plots[[1L]])
    if (isTRUE(combined)) return(.combine_plot_grid(plots, fallback = "input"))
    plots
}

.expand_association_terms <- function(requested_terms, available_terms) {
    # Normalize both inputs early so later matching works on one canonical
    # character representation regardless of how the caller supplied terms.
    requested_terms <- unique(as.character(requested_terms %||% character(0)))
    available_terms <- unique(as.character(available_terms %||% character(0)))
    if (!length(requested_terms) || !length(available_terms)) {
        return(character(0))
    }

    # Expand a channel request like "vcov" to every concrete component term
    # such as vcov[1], vcov[2], ... while still allowing exact one-to-one
    # requests for already-expanded names.
    matched <- unlist(lapply(requested_terms, function(requested_term) {
        hits <- available_terms[
            available_terms == requested_term |
                startsWith(available_terms, paste0(requested_term, "["))
        ]
        hits
    }), use.names = FALSE)

    unique(matched)
}

.available_association_terms <- function(x) {
    # The cached data is the preferred source because it already knows the
    # exact user-facing term labels after any corr/vcov component expansion.
    data <- .get_association_plot_data(x)
    if (!is.null(data$term_map) && nrow(data$term_map) > 0L) {
        terms <- unique(as.character(data$term_map$term))
        return(terms[!grepl("^weight:\\s", terms)])
    }

    # Fall back to a fresh association-draw extraction only when the data is
    # absent; this keeps plotting robust even if caching was skipped earlier.
    assoc_draws <- extract.JoinMeFit(x, what = "assoc", keep_chains = FALSE)
    terms <- unique(as.character(assoc_draws$term_map$term))
    terms[!grepl("^weight:\\s", terms)]
}

.plot_joinmefit_association_single <- function(x, term, marker, association_grid, association_range, association_points, ci_levels,
                                               ci_type, association_metric, prediction_style, theme_fn, seed) {
    # Map concrete component labels back to the underlying transform channel so
    # downstream helpers can reuse corr/vcov-specific machinery while still
    # preserving the exact plotted component term.
    term_key <- if (grepl("^corr", term)) {
        "corr"
    } else if (grepl("^vcov", term)) {
        "vcov"
    } else {
        term
    }

    # Reuse one shared data bundle so every helper in this plotting pass sees the
    # same cached support, term map, transform draws, and marker weights.
    data <- .get_association_plot_data(x, seed = seed)

    # Pull the association-coefficient draws for exactly the requested plotted
    # term. For component-expanded corr/vcov terms this returns the single
    # column matching that component, not the whole channel matrix.
    coeff_draws <- .assoc_coeff_draws(x, term = term, data = data, seed = seed)

    # Determine whether the term should be plotted once globally (corr/vcov) or
    # once per marker-level trajectory-style channel.
    marker_levels <- .association_plot_markers(x, term_key, marker)

    # Precompute quantile probabilities and output column names once so the
    # per-marker loop can focus only on generating posterior curves.
    probs <- .quantile_probs_from_ci_plot(ci_levels)
    q_names <- .quantile_colnames(probs)

    # Build one plotting data block per marker, then stack them. Each block is
    # generated from the same algorithm: choose x grid -> transform x grid for
    # every draw -> optionally reweight -> optionally convert to hazard-scale
    # contribution -> summarize draw-wise uncertainty at each x location.
    plot_df <- do.call(rbind, lapply(marker_levels, function(marker_level) {
        # Either respect an explicit user grid or synthesize one from model-
        # implied support for the exact plotted term/component.
        x_grid <- association_grid %||% .default_association_grid(
            x,
            term_key,
            term = term,
            marker = marker_level,
            association_range = association_range,
            n = association_points
        )

        # Force a clean sorted numeric grid because later quantile summaries and
        # ggplot layers assume monotone x values with duplicates removed.
        x_grid <- sort(unique(as.numeric(x_grid)))

        # Evaluate the association transform on the raw x-grid for every draw.
        # The returned matrix has one row per posterior draw and one column per
        # x location, which keeps later hazard-vs-transform branching simple.
        tf_mat <- .association_transform_matrix(x, term_key, term = term, x_grid = x_grid, n_draws = length(coeff_draws), seed = seed, data = data)

        # Marker-weighted channels such as cv_total/cs_total need one extra
        # multiplicative layer so the plotted transform matches the weighted
        # marker aggregation used during fitting and prediction.
        weight_draws <- .association_marker_weight_draws(
            x = x,
            term_key = term_key,
            marker_level = marker_level,
            n_draws = length(coeff_draws),
            seed = seed,
            data = data
        )
        tf_mat <- tf_mat * as.numeric(weight_draws)

        # Hazard-scale corr/vcov contributions are zero-referenced at raw value
        # 0 so the plotted curve shows only the incremental contribution beyond
        # the baseline hazard, matching the fitted model semantics.
        if (identical(association_metric, "hazard") && term_key %in% c("corr", "vcov")) {
            tf_ref <- .association_transform_matrix(
                x,
                term_key,
                term = term,
                x_grid = 0,
                n_draws = length(coeff_draws),
                seed = seed,
                data = data
            )
            tf_mat <- tf_mat - matrix(tf_ref[, 1], nrow = nrow(tf_mat), ncol = ncol(tf_mat))
        }

        # A transform plot shows f(x) directly, while a hazard plot multiplies
        # the transformed values by the posterior association coefficient draws.
        curve_mat <- if (identical(association_metric, "hazard")) {
            tf_mat * coeff_draws
        } else {
            tf_mat
        }

        # Summarize the posterior curve column-by-column so each x position gets
        # mean, sd, and all requested quantile bands.
        q_mat <- t(apply(curve_mat, 2, stats::quantile, probs = probs, na.rm = TRUE, names = FALSE))
        out <- data.frame(
            x = x_grid,
            mean = colMeans(curve_mat),
            sd = apply(curve_mat, 2, stats::sd),
            marker = as.character(marker_level),
            stringsAsFactors = FALSE
        )
        out[q_names] <- q_mat
        out
    }))

    # Use the posterior median as the main central trajectory, falling back to
    # the already-computed quantile naming convention used elsewhere in plotting.
    median_col <- .quantile_name_from_prob(0.5)

    # A single corr/vcov component is plotted globally, but mean/current-value
    # channels may produce one curve per marker; the aesthetics change slightly
    # depending on which case applies.
    multi_marker <- length(unique(plot_df$marker)) > 1L
    p <- if (multi_marker) {
        ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data[[median_col]], color = .data$marker, fill = .data$marker, group = .data$marker))
    } else {
        ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data[[median_col]]))
    }

    # Add ribbons from widest interval to narrowest so narrower credible bands
    # are drawn on top and remain visible.
    if (ci_type %in% c("ribbon", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            nm <- .quantile_names_from_ci(level)
            if (all(nm %in% names(plot_df))) {
                p <- if (multi_marker) {
                    p + ggplot2::geom_ribbon(
                        ggplot2::aes(ymin = .data[[nm[1]]], ymax = .data[[nm[2]]]),
                        alpha = (prediction_style$alpha %||% 0.2) / length(ci_levels),
                        color = NA
                    )
                } else {
                    p + ggplot2::geom_ribbon(
                        ggplot2::aes(ymin = .data[[nm[1]]], ymax = .data[[nm[2]]]),
                        fill = prediction_style$fill %||% "steelblue",
                        alpha = (prediction_style$alpha %||% 0.2) / length(ci_levels),
                        color = NA
                    )
                }
            }
        }
    }

    # Draw the central curve after ribbons so it remains visually prominent.
    p <- if (multi_marker) {
        p + ggplot2::geom_line(linewidth = prediction_style$linewidth %||% 0.8)
    } else {
        p + ggplot2::geom_line(
            color = prediction_style$color %||% "steelblue",
            linewidth = prediction_style$linewidth %||% 0.8
        )
    }

    # Optional dashed quantile boundaries provide an alternative CI style or a
    # supplement to the ribbons when `ci_type = "both"`.
    if (ci_type %in% c("line", "both")) {
        for (level in sort(ci_levels, decreasing = TRUE)) {
            nm <- .quantile_names_from_ci(level)
            for (qcol in nm) {
                if (qcol %in% names(plot_df)) {
                    p <- if (multi_marker) {
                        p + ggplot2::geom_line(
                            ggplot2::aes(y = .data[[qcol]]),
                            linetype = "dashed",
                            alpha = 0.5,
                            linewidth = 0.4
                        )
                    } else {
                        p + ggplot2::geom_line(
                            ggplot2::aes(y = .data[[qcol]]),
                            color = prediction_style$color %||% "steelblue",
                            linetype = "dashed",
                            alpha = 0.5,
                            linewidth = 0.4
                        )
                    }
                }
            }
        }
    }

    # Use different labels for hazard-scale contributions vs pure transforms so
    # the viewer can immediately tell whether the plot is on the hazard scale
    # or shows the transform alone.
    y_lab <- if (identical(association_metric, "hazard")) {
        latex2exp::TeX("$\\delta \\log H$")
    } else {
        paste0("f(", term, ")")
    }

    plot_title <- if (identical(association_metric, "hazard")) {
        paste("Association contribution:", term)
    } else {
        paste("Association transform:", term)
    }

    # Finalize plot metadata after all statistical content is assembled.
    p <- p + ggplot2::labs(
        title = plot_title,
        x = term,
        y = y_lab,
        color = if (multi_marker) "marker" else NULL,
        fill = if (multi_marker) "marker" else NULL
    ) + theme_fn()

    # Some transform-only channels naturally live on a bounded support like
    # [0, 1]; apply the corresponding y-scale adjustment only at the end so it
    # sees the full summarized plotting data.
    y_scale_adjust <- .assoc_transform_limits(plot_df, association_metric)
    if (!is.null(y_scale_adjust)) {
        p <- p + y_scale_adjust
    }
    p
}

.association_plot_markers <- function(x, term_key, marker = NULL) {
    # corr/vcov channels are global subject-level features, so marker faceting
    # would be misleading; force a single synthetic marker label in that case.
    if (term_key %in% c("corr", "vcov")) {
        return("all")
    }

    # For marker-resolved channels, discover the available marker levels from
    # the fitted longitudinal data so the plot mirrors the original fit input.
    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    available_markers <- if (!is.null(x$dataLong) && marker_var %in% names(x$dataLong)) {
        unique(as.character(stats::na.omit(x$dataLong[[marker_var]])))
    } else {
        character(0)
    }

    # If the user did not request a subset, plot all available markers; if the
    # data carry no marker labels, fall back to one synthetic "all" marker.
    if (is.null(marker)) {
        if (length(available_markers) == 0) "all" else available_markers
    } else {
        marker
    }
}

.transform_specs <- function(x) {
    # Keep transform lookup centralized because some fits store the resolved
    # transform spec in config while others only retain the original call.
    x$config$transforms_spec %||% x$call$transforms %||% list()
}

.assoc_component_index <- function(term) {
    # Most channels are scalar and therefore default to component 1. Expanded
    # corr/vcov labels encode the component index in square brackets.
    term <- as.character(term %||% "")[1]
    if (!grepl("\\[\\d+\\]$", term)) {
        return(1L)
    }
    as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", term))
}

.transform_spec_for_term <- function(x, term_key, term = term_key) {
    # The current plotting API shares one transform spec per channel, even when
    # corr/vcov expand into multiple component terms.
    specs <- .transform_specs(x)
    specs[[term_key]] %||% list(type = "identity")
}

.plot_raw_knot_range <- function(tf_spec) {
    # Recover the raw plotting range implied by either explicit knots or pwlin
    # x-values. This is used as a safe fallback whenever observed/model support
    # reaches outside the fitted spline domain.
    knots <- as.numeric(tf_spec$knots %||% tf_spec$x %||% numeric(0))
    if (!length(knots) || !all(is.finite(knots))) {
        return(NULL)
    }

    # Expit-based spline transforms store knot locations on the bounded expit
    # domain, but plotting works on the raw pre-transform scale, so convert back
    # with qlogis after guarding against exact 0/1 endpoints.
    if (.transform_uses_expit_input(tf_spec)) {
        knots <- .validate_expit_domain_values(knots, "knots")
        eps <- sqrt(.Machine$double.eps)
        knots <- stats::qlogis(pmin(pmax(knots, eps), 1 - eps))
    }
    if (!all(is.finite(knots))) {
        return(NULL)
    }
    range(knots)
}

.support_outside_spline_range <- function(support, tf_spec) {
    # Determine whether a proposed raw support interval would require spline
    # extrapolation beyond the fitted boundary knots.
    knots <- as.numeric(tf_spec$knots %||% numeric(0))
    if (length(support) != 2L || length(knots) < 2L) {
        return(FALSE)
    }
    boundary <- c(knots[1], knots[length(knots)])

    # Expit-input transforms compare support on the bounded expit scale, not on
    # the raw scale used elsewhere in plotting.
    if (.transform_uses_expit_input(tf_spec)) {
        support <- .transform_input_for_spec(support, tf_spec)
    }
    any(support < boundary[1] | support > boundary[2])
}

.default_association_grid <- function(x, term_key, term = term_key, marker = NULL, association_range = NULL, n = 200) {
    # An explicit user range always wins; all other heuristics are only for the
    # default auto-grid construction.
    if (!is.null(association_range)) {
        return(seq(association_range[1], association_range[2], length.out = n))
    }

    # corr has a fixed theoretical raw support of [-1, 1], but if model-implied
    # component-specific support is cached we prefer that narrower interval.
    if (identical(term_key, "corr")) {
        support <- .association_support_range(x, term_key = term_key, term = term, marker = marker)
        if (!is.null(support)) {
            return(seq(support[1], support[2], length.out = n))
        }
        return(seq(-1, 1, length.out = n))
    }

    # vcov does not have a universal theoretical support because it mixes SD and
    # off-diagonal Cholesky features, so rely on cached model-implied support.
    if (identical(term_key, "vcov")) {
        support <- .association_support_range(x, term_key = term_key, term = term, marker = marker)
        if (!is.null(support)) {
            return(seq(support[1], support[2], length.out = n))
        }
    }

    # For all other channels, consult the transform spec so auto-grid selection
    # can avoid invalid spline extrapolation when support and knot ranges clash.
    tf_spec <- .transform_specs(x)[[term_key]]
    tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")
    support <- .association_support_range(x, term_key = term_key, term = term, marker = marker)
    if (!is.null(support)) {
        if (.is_ispline_transform_type(tf_type) && !is.null(tf_spec$knots)) {
            kr <- .plot_raw_knot_range(tf_spec)
            if (!is.null(kr) && diff(kr) > 0 && .support_outside_spline_range(support, tf_spec)) {
                cli::cli_warn(c(
                    x = "Model-implied support for {.val {term_key}} extends beyond the fitted spline knot range.",
                    i = "Using knot support [{format(signif(kr[1], 4), scientific = FALSE)}, {format(signif(kr[2], 4), scientific = FALSE)}] to avoid unsupported spline extrapolation."
                ))
                return(seq(kr[1], kr[2], length.out = n))
            }
        }
        return(seq(support[1], support[2], length.out = n))
    }

    # If there is no cached support, fall back to heuristics based on observed
    # data, starting from the relevant marker subset when applicable.
    response_var <- all.vars(x$formulaLong)[1] %||% "y"
    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    data_long <- x$dataLong
    if (!is.null(marker) && marker_var %in% names(data_long)) {
        data_long <- data_long[as.character(data_long[[marker_var]]) %in% as.character(marker), , drop = FALSE]
    }
    y_obs <- data_long[[response_var]]

    # Slope channels derive their natural raw support from observed finite-
    # difference slopes within each (id, marker) trajectory.
    if (term_key %in% c("cs_total", "cs_mean", "cs_marker")) {
        id_var <- .joinmefit_call_arg_chr(x$call, "id_var", "id")
        time_var <- .joinmefit_call_arg_chr(x$call, "time_var", "time")
        dat <- data_long[, c(id_var, time_var, marker_var, response_var), drop = FALSE]
        names(dat) <- c("id", "time", "marker", "y")
        dat <- dat[order(dat$id, dat$marker, dat$time), , drop = FALSE]
        dy <- ave(dat$y, interaction(dat$id, dat$marker), FUN = function(v) c(NA_real_, diff(v)))
        dt <- ave(dat$time, interaction(dat$id, dat$marker), FUN = function(v) c(NA_real_, diff(v)))
        slope <- dy / dt
        slope <- slope[is.finite(slope)]
        if (length(slope) > 1L) {
            xr <- stats::quantile(slope, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
            if (all(is.finite(xr)) && diff(xr) > 0) return(seq(xr[1], xr[2], length.out = n))
        }
    }

    # Current-value channels use observed response quantiles as a pragmatic raw
    # support estimate unless spline-knot safety forces a narrower range.
    if (term_key %in% c("cv_total", "cv_mean", "cv_marker") && length(y_obs) > 1L) {
        xr <- stats::quantile(y_obs, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
        if (all(is.finite(xr)) && diff(xr) > 0) {
            if (.is_ispline_transform_type(tf_type) && !is.null(tf_spec$knots)) {
                kr <- .plot_raw_knot_range(tf_spec)
                if (!is.null(kr) && diff(kr) > 0 && .support_outside_spline_range(xr, tf_spec)) {
                    cli::cli_warn(c(
                        x = "Observed support for {.val {term_key}} extends beyond the fitted spline knot range.",
                        i = "Using knot support [{format(signif(kr[1], 4), scientific = FALSE)}, {format(signif(kr[2], 4), scientific = FALSE)}] to avoid unsupported spline extrapolation."
                    ))
                    return(seq(kr[1], kr[2], length.out = n))
                }
            }
            return(seq(xr[1], xr[2], length.out = n))
        }
    }

    # If user-specified transform inputs define a natural raw domain, prefer
    # that over the generic final fallback.
    if (!is.null(tf_spec$x)) {
        xr <- .plot_raw_knot_range(list(type = tf_type, x = tf_spec$x))
        if (all(is.finite(xr)) && diff(xr) > 0) {
            return(seq(xr[1], xr[2], length.out = n))
        }
    }
    if (!is.null(tf_spec$knots)) {
        xr <- .plot_raw_knot_range(tf_spec)
        if (all(is.finite(xr)) && diff(xr) > 0) {
            return(seq(xr[1], xr[2], length.out = n))
        }
    }

    # As a last empirical fallback, reuse observed response quantiles even for
    # channels without a more specialized support heuristic.
    if (length(y_obs) > 1L) {
        xr <- stats::quantile(y_obs, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
        if (all(is.finite(xr)) && diff(xr) > 0) return(seq(xr[1], xr[2], length.out = n))
    }

    # The generic symmetric range only applies when every data- and model-based
    # heuristic above failed to identify a more informative support interval.
    seq(-1, 1, length.out = n)
}

.assoc_transform_limits <- function(plot_df, association_metric) {
    # Only transform plots use this bounded-scale heuristic; hazard plots should
    # retain their natural contribution scale.
    if (!identical(association_metric, "transform")) {
        return(NULL)
    }

    # Inspect the full set of summarized y-values, not just the median, so the
    # decision reflects all visible uncertainty bands.
    q_cols <- grep("^q", names(plot_df), value = TRUE)
    vals <- unlist(plot_df[, unique(c("mean", q_cols)), drop = FALSE], use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
        return(NULL)
    }

    # When the transform appears to live on an approximate probability scale,
    # clamp the displayed y-axis accordingly for readability.
    if (min(vals) >= -0.02 && max(vals) <= 1.02) {
        return(ggplot2::coord_cartesian(ylim = c(0, 1)))
    }

    NULL
}

.association_marker_weight_draws <- function(x, term_key, marker_level, n_draws, seed = 1, data = NULL) {
    # Most association channels are not marker-weighted, so return an identity
    # multiplier and let the caller reuse one code path for all channels.
    if (!(term_key %in% c("cv_total", "cs_total", "cv_marker", "cs_marker"))) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    # Discover the marker ordering used by standata so weight draws line up with
    # the same marker index convention as the fitted model.
    marker_levels <- as.character(x$stan_data$marker_levels %||% unique(stats::na.omit(x$dataLong$marker)) %||% character(0))
    if (length(marker_levels) == 0L) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    marker_idx <- match(as.character(marker_level), marker_levels)
    if (is.na(marker_idx)) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    weight_draws <- .marker_weight_draws(x, term_key = term_key, n_draws = n_draws, seed = seed, data = data)
    if (is.null(weight_draws) || ncol(weight_draws) < marker_idx) {
        return(matrix(1 / length(marker_levels), nrow = n_draws, ncol = 1L))
    }

    # Divide by the number of markers because the cv/cs total-style channels are
    # plotted on the same weighted-average scale used elsewhere in the package.
    matrix(weight_draws[, marker_idx, drop = TRUE] / length(marker_levels), ncol = 1L)
}

.marker_weight_draws <- function(x, term_key, n_draws, seed = 1, data = NULL) {
    # Start from the canonical marker ordering used in standata so cached draws,
    # posterior draws, and default weights all align to the same columns.
    sd <- x$stan_data
    marker_levels <- as.character(sd$marker_levels %||% unique(stats::na.omit(x$dataLong$marker)) %||% character(0))
    n_markers <- length(marker_levels)
    if (n_markers == 0L) {
        return(NULL)
    }

    # Prefer cached data draws because they may already be subset to the
    # relevant variables and draws for the current plotting request.
    if (!is.null(data$marker_weight_draws)) {
        if (is.matrix(data$marker_weight_draws)) {
            return(as.matrix(data$marker_weight_draws[, seq_len(min(n_markers, ncol(data$marker_weight_draws))), drop = FALSE]))
        }
        cached <- data$marker_weight_draws[[term_key]]
        if (!is.null(cached)) {
            return(as.matrix(cached[, seq_len(min(n_markers, ncol(cached))), drop = FALSE]))
        }
    }

    # Otherwise try the effective-weight parameter names first.
    eff_names <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(n_markers), "]")
    dmat <- tryCatch(
        .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
        error = function(e) NULL
    )
    if ((is.null(dmat) || !all(eff_names %in% colnames(dmat))) && isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)) {
        eff_names <- paste0("marker_weights_eff[", seq_len(n_markers), "]")
        dmat <- tryCatch(
            .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
            error = function(e) NULL
        )
    }
    if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
        return(as.matrix(dmat[, eff_names, drop = FALSE]))
    }

    # Fall back to the pre-effective base-weight parameterization if that is all
    # that is available in the stored fit object.
    base_names <- paste0(.marker_weight_var_prefix(term_key, effective = FALSE), "[", seq_len(n_markers), "]")
    dmat <- tryCatch(
        .get_draws_matrix(x$fit, variables = base_names, draws = n_draws, seed = seed),
        error = function(e) NULL
    )
    if ((is.null(dmat) || !all(base_names %in% colnames(dmat))) && isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)) {
        base_names <- paste0("marker_weights[", seq_len(n_markers), "]")
        dmat <- tryCatch(
            .get_draws_matrix(x$fit, variables = base_names, draws = n_draws, seed = seed),
            error = function(e) NULL
        )
    }
    if (!is.null(dmat) && all(base_names %in% colnames(dmat))) {
        return(as.matrix(dmat[, base_names, drop = FALSE]))
    }

    # If no posterior draws are available, degrade gracefully to the standata
    # base weights so plotting still works in lightweight or partial objects.
    base_weights <- as.numeric((sd$marker_weights_by_term %||% list())[[term_key]] %||% sd$marker_weights %||% rep(1, n_markers))
    if (length(base_weights) < n_markers) {
        base_weights <- c(base_weights, rep(1, n_markers - length(base_weights)))
    }
    matrix(rep(base_weights[seq_len(n_markers)], each = n_draws), nrow = n_draws, byrow = FALSE)
}

.assoc_channel_map <- function(term_key) {
    # Centralize the mapping from public association channel names to the Stan
    # variable prefixes and spline metadata used when reconstructing transform
    # matrices from posterior draws.
    switch(term_key,
        cv_total = list(
            eff_prefix = "coeff_cv_eff",
            base_coeff = "coeff_cv",
            n_coeff = "n_coeff_cv",
            knots = "knots_cv",
            degree = "spline_degree_cv",
            iota_intercept = "iota_intercept_cv_eff",
            iota_slope = "iota_slope_cv_eff",
            n_iota_intercept = "estimate_iota_intercept_cv",
            n_iota_slope = "estimate_iota_slope_cv"
        ),
        cs_total = list(
            eff_prefix = "coeff_cs_eff",
            base_coeff = "coeff_cs",
            n_coeff = "n_coeff_cs",
            knots = "knots_cs",
            degree = "spline_degree_cs",
            iota_intercept = "iota_intercept_cs_eff",
            iota_slope = "iota_slope_cs_eff",
            n_iota_intercept = "estimate_iota_intercept_cs",
            n_iota_slope = "estimate_iota_slope_cs"
        ),
        corr = list(
            eff_prefix = "coeff_corr_eff",
            base_coeff = "coeff_corr",
            n_coeff = "n_coeff_corr",
            knots = "knots_corr",
            degree = "spline_degree_corr",
            iota_intercept = "iota_intercept_corr_eff",
            iota_slope = "iota_slope_corr_eff",
            n_iota_intercept = "estimate_iota_intercept_corr",
            n_iota_slope = "estimate_iota_slope_corr"
        ),
        vcov = list(
            eff_prefix = "coeff_vcov_eff",
            base_coeff = "coeff_vcov",
            n_coeff = "n_coeff_vcov",
            knots = "knots_vcov",
            degree = "spline_degree_vcov",
            iota_intercept = "iota_intercept_vcov_eff",
            iota_slope = "iota_slope_vcov_eff",
            n_iota_intercept = "estimate_iota_intercept_vcov",
            n_iota_slope = "estimate_iota_slope_vcov"
        ),
        cv_mean = list(
            eff_prefix = "coeff_cv_mean_eff",
            base_coeff = "coeff_cv_mean",
            n_coeff = "n_coeff_cv_mean",
            knots = "knots_cv_mean",
            degree = "spline_degree_cv_mean",
            iota_intercept = "iota_intercept_cv_mean_eff",
            iota_slope = "iota_slope_cv_mean_eff",
            n_iota_intercept = "estimate_iota_intercept_cv_mean",
            n_iota_slope = "estimate_iota_slope_cv_mean"
        ),
        cv_marker = list(
            eff_prefix = "coeff_cv_marker_eff",
            base_coeff = "coeff_cv_marker",
            n_coeff = "n_coeff_cv_marker",
            knots = "knots_cv_marker",
            degree = "spline_degree_cv_marker",
            iota_intercept = "iota_intercept_cv_marker_eff",
            iota_slope = "iota_slope_cv_marker_eff",
            n_iota_intercept = "estimate_iota_intercept_cv_marker",
            n_iota_slope = "estimate_iota_slope_cv_marker"
        ),
        cs_mean = list(
            eff_prefix = "coeff_cs_mean_eff",
            base_coeff = "coeff_cs_mean",
            n_coeff = "n_coeff_cs_mean",
            knots = "knots_cs_mean",
            degree = "spline_degree_cs_mean",
            iota_intercept = "iota_intercept_cs_mean_eff",
            iota_slope = "iota_slope_cs_mean_eff",
            n_iota_intercept = "estimate_iota_intercept_cs_mean",
            n_iota_slope = "estimate_iota_slope_cs_mean"
        ),
        cs_marker = list(
            eff_prefix = "coeff_cs_marker_eff",
            base_coeff = "coeff_cs_marker",
            n_coeff = "n_coeff_cs_marker",
            knots = "knots_cs_marker",
            degree = "spline_degree_cs_marker",
            iota_intercept = "iota_intercept_cs_marker_eff",
            iota_slope = "iota_slope_cs_marker_eff",
            n_iota_intercept = "estimate_iota_intercept_cs_marker",
            n_iota_slope = "estimate_iota_slope_cs_marker"
        ),
        NULL
    )
}

.association_transform_matrix <- function(x, term_key, term = term_key, x_grid, n_draws, seed = 1, data = NULL) {
    # Resolve the transform specification for the requested channel/component,
    # then route to the lightest-weight evaluator for that transform family.
    tf_spec <- .transform_spec_for_term(x, term_key, term = term)
    tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")

    # Identity transforms are the cheapest case: replicate the raw grid across
    # posterior draws without any additional computation.
    if (identical(tf_type, "identity")) {
        return(matrix(rep(as.numeric(x_grid), each = n_draws), nrow = n_draws))
    }

    # Functional and pwlin transforms can be evaluated once on the grid because
    # they do not depend on draw-specific spline coefficients.
    if (identical(tf_type, "functional")) {
        tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
        iota_draws <- .transform_iota_draws(
            x = x,
            term_key = term_key,
            term = term,
            n_draws = n_draws,
            seed = seed,
            data = data
        )
        intercept_draws <- as.matrix(iota_draws$intercept %||% matrix(0, nrow = n_draws, ncol = 0L))
        slope_draws <- as.matrix(iota_draws$slope %||% matrix(1, nrow = n_draws, ncol = 0L))
        if (nrow(intercept_draws) != n_draws) {
            intercept_draws <- matrix(rep(intercept_draws[1, , drop = TRUE], each = n_draws), nrow = n_draws)
        }
        if (nrow(slope_draws) != n_draws) {
            slope_draws <- matrix(rep(slope_draws[1, , drop = TRUE], each = n_draws), nrow = n_draws)
        }
        shifted <- vapply(seq_len(n_draws), function(i) {
            as.numeric(tf_fun(
                x_grid,
                iota_intercept = intercept_draws[i, , drop = TRUE],
                iota_slope = slope_draws[i, , drop = TRUE]
            ))
        }, numeric(length(x_grid)))
        return(t(shifted))
    }

    if (identical(tf_type, "pwlin")) {
        tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
        vals <- as.numeric(tf_fun(x_grid))
        return(matrix(rep(vals, each = n_draws), nrow = n_draws))
    }

    # Spline transforms do depend on posterior coefficient draws, so they need
    # their dedicated matrix builder.
    if (.is_ispline_transform_type(tf_type)) {
        return(.ispline_transform_matrix(x, term_key, term = term, x_grid, n_draws = n_draws, seed = seed, data = data))
    }

    # The final fallback keeps plotting resilient to future transform types that
    # can still be represented by a deterministic pointwise transform function.
    tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
    vals <- as.numeric(tf_fun(x_grid))
    matrix(rep(vals, each = n_draws), nrow = n_draws)
}

.ispline_transform_matrix <- function(x, term_key, term = term_key, x_grid, n_draws, seed = 1, data = NULL) {
    if (!requireNamespace("splines2", quietly = TRUE)) {
        cli::cli_abort(c(
            x = "Package {.pkg splines2} is required for spline-based association plotting.",
            i = "Install {.pkg splines2} to plot spline transforms."
        ))
    }

    sd <- x$stan_data
    map <- .assoc_channel_map(term_key)
    tf_spec <- .transform_spec_for_term(x, term_key, term = term)

    knots <- as.numeric(sd[[map$knots]] %||% tf_spec$knots %||% tf_spec$x)
    degree <- as.integer(sd[[map$degree]] %||% tf_spec$degree %||% 3L)
    n_coeff <- as.integer(sd[[map$n_coeff]] %||% length(sd[[map$base_coeff]] %||% tf_spec$coeff %||% numeric(0)))
    coeff_draws <- .transform_coeff_draws(
        x = x,
        eff_prefix = map$eff_prefix,
        base_coeff = sd[[map$base_coeff]] %||% tf_spec$coeff,
        n_coeff = n_coeff,
        component_index = if (term_key %in% c("corr", "vcov")) .assoc_component_index(term) else NULL,
        n_draws = n_draws,
        seed = seed,
        data = data,
        term_key = term_key,
        term = term
    )

    if (length(knots) < 2L) {
        cli::cli_abort(c(
            x = "Spline transform for {.val {term_key}} is missing knot information.",
            i = "Refit with stored spline metadata or provide an explicit transform specification."
        ))
    }

    boundary <- c(knots[1], knots[length(knots)])
    internal_knots <- if (length(knots) > 2L) knots[2:(length(knots) - 1L)] else numeric(0)
    x_eval_basis <- .transform_input_for_spec(as.numeric(x_grid), tf_spec)
    x_eval_basis <- pmin(pmax(x_eval_basis, boundary[1]), boundary[2])
    basis <- splines2::iSpline(
        x_eval_basis,
        knots = internal_knots,
        degree = degree,
        intercept = TRUE,
        Boundary.knots = boundary
    )
    basis <- as.matrix(basis)

    n_basis <- ncol(basis)
    if (ncol(coeff_draws) < n_basis) {
        coeff_draws <- cbind(coeff_draws, matrix(0, nrow(coeff_draws), n_basis - ncol(coeff_draws)))
    } else if (ncol(coeff_draws) > n_basis) {
        coeff_draws <- coeff_draws[, seq_len(n_basis), drop = FALSE]
    }

    out <- coeff_draws %*% t(basis)
    if (.resolve_monotone_direction(tf_spec$direction %||% tf_spec$spline_direction %||% 1L, default = 1L) < 0L) {
        out <- matrix(rowSums(coeff_draws), nrow = nrow(coeff_draws), ncol = ncol(out)) - out
    }
    out
}

.transform_coeff_draws <- function(x, eff_prefix, base_coeff, n_coeff, n_draws, seed = 1, data = NULL, term_key = NULL, term = term_key, component_index = NULL) {
    n_coeff <- as.integer(n_coeff %||% 0L)
    if (n_coeff < 1L) {
        return(matrix(0, nrow = n_draws, ncol = 0L))
    }

    data_key <- term %||% term_key
    if (!is.null(data_key) && !is.null(data$transform_coeff_draws[[data_key]])) {
        return(as.matrix(data$transform_coeff_draws[[data_key]]))
    }

    if (!is.null(component_index) && !is.na(component_index)) {
        eff_names <- paste0(eff_prefix, "[", component_index, ",", seq_len(n_coeff), "]")
        dmat <- tryCatch(
            .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
            error = function(e) NULL
        )
        if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
            return(as.matrix(dmat[, eff_names, drop = FALSE]))
        }

        base_mat <- as.matrix(base_coeff %||% matrix(0, nrow = component_index, ncol = n_coeff))
        if (nrow(base_mat) < component_index) {
            base_mat <- rbind(base_mat, matrix(0, nrow = component_index - nrow(base_mat), ncol = ncol(base_mat)))
        }
        if (ncol(base_mat) < n_coeff) {
            base_mat <- cbind(base_mat, matrix(0, nrow = nrow(base_mat), ncol = n_coeff - ncol(base_mat)))
        }
        base_vec <- base_mat[component_index, seq_len(n_coeff), drop = TRUE]
        return(matrix(rep(base_vec, times = n_draws), nrow = n_draws, byrow = TRUE))
    }

    eff_names <- paste0(eff_prefix, "[", seq_len(n_coeff), "]")
    dmat <- tryCatch(
        .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
        error = function(e) NULL
    )
    if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
        return(as.matrix(dmat[, eff_names, drop = FALSE]))
    }

    base_coeff <- as.numeric(base_coeff %||% rep(0, n_coeff))
    if (length(base_coeff) < n_coeff) {
        base_coeff <- c(base_coeff, rep(0, n_coeff - length(base_coeff)))
    }
    matrix(rep(base_coeff[seq_len(n_coeff)], times = n_draws), nrow = n_draws, byrow = TRUE)
}

.transform_iota_draws <- function(x, term_key, term = term_key, n_draws, seed = 1, data = NULL) {
    data_key <- term %||% term_key
    if (!is.null(data_key) && !is.null(data$transform_iota_draws[[data_key]])) {
        cached <- data$transform_iota_draws[[data_key]]
        return(list(
            intercept = as.matrix(cached$intercept %||% matrix(0, nrow = n_draws, ncol = 0L)),
            slope = as.matrix(cached$slope %||% matrix(1, nrow = n_draws, ncol = 0L))
        ))
    }

    tf_spec <- .transform_spec_for_term(x, term_key, term = term)
    if (!identical(.canonicalise_transform_type(tf_spec$type %||% "identity"), "functional")) {
        return(list(
            intercept = matrix(0, nrow = n_draws, ncol = 0L),
            slope = matrix(1, nrow = n_draws, ncol = 0L)
        ))
    }

    map <- .assoc_channel_map(term_key)
    component_index <- if (term_key %in% c("corr", "vcov")) .assoc_component_index(term) else NULL
    n_iota_intercept <- as.integer(tf_spec$n_iota_intercept %||% x$stan_data[[map$n_iota_intercept]] %||% 0L)
    n_iota_slope <- as.integer(tf_spec$n_iota_slope %||% x$stan_data[[map$n_iota_slope]] %||% 0L)

    get_iota_component <- function(prefix, default, n_iota) {
        if (is.null(prefix) || n_iota < 1L) {
            return(matrix(default, nrow = n_draws, ncol = max(0L, n_iota)))
        }
        if (!is.null(component_index) && !is.na(component_index)) {
            offsets <- (component_index - 1L) * n_iota + seq_len(n_iota)
            eff_names <- paste0(prefix, "[", offsets, "]")
        } else {
            eff_names <- paste0(prefix, "[", seq_len(n_iota), "]")
        }
        dmat <- tryCatch(
            .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
            error = function(e) NULL
        )
        if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
            return(as.matrix(dmat[, eff_names, drop = FALSE]))
        }
        if (n_iota == 1L && is.null(component_index) && prefix %in% posterior::variables(.get_draws_obj(x$fit))) {
            dmat_single <- tryCatch(
                .get_draws_matrix(x$fit, variables = prefix, draws = n_draws, seed = seed),
                error = function(e) NULL
            )
            if (!is.null(dmat_single) && prefix %in% colnames(dmat_single)) {
                return(matrix(as.numeric(dmat_single[, prefix, drop = TRUE]), ncol = 1L))
            }
        }
        matrix(default, nrow = n_draws, ncol = n_iota)
    }

    list(
        intercept = get_iota_component(map$iota_intercept, default = 0, n_iota = n_iota_intercept),
        slope = get_iota_component(map$iota_slope, default = 1, n_iota = n_iota_slope)
    )
}

.assoc_coeff_draws <- function(x, term, data = NULL, seed = 1) {
    term_key <- if (grepl("^corr", term)) {
        "corr"
    } else if (grepl("^vcov", term)) {
        "vcov"
    } else {
        term
    }
    if (!is.null(data$coeff_draws[[term_key]])) {
        vals <- data$coeff_draws[[term_key]]
        if (is.matrix(vals)) {
            target_var_candidates <- character(0)
            if (!is.null(data$term_map) && nrow(data$term_map) > 0L) {
                target_rows <- data$term_map[data$term_map$term == term, , drop = FALSE]
                if (nrow(target_rows) > 0L) {
                    target_var_candidates <- as.character(target_rows$variable)
                }
            }
            if (grepl("^(corr|vcov)\\[\\d+\\]$", term)) {
                target_var_candidates <- unique(c(
                    target_var_candidates,
                    paste0("alpha_", sub("\\[", "_eff[", term)),
                    paste0("alpha_", term)
                ))
            }
            if (length(target_var_candidates) > 0L && !is.null(colnames(vals))) {
                matched_var <- target_var_candidates[target_var_candidates %in% colnames(vals)][1L]
                if (!is.na(matched_var) && nzchar(matched_var)) {
                    return(as.numeric(vals[, matched_var, drop = TRUE]))
                }
            }
            return(as.numeric(vals[, 1]))
        }
        return(as.numeric(vals))
    }

    assoc_draws <- extract.JoinMeFit(x, what = "assoc", term = term, keep_chains = FALSE)$draws
    as.numeric(assoc_draws[, 1])
}

.get_association_plot_data <- function(x, seed = 1) {
    data <- x$config$association_plot_data %||% x$cache_get("association_plot_data")
    if (!is.null(data)) {
        return(data)
    }
    if (is.null(x$fit)) {
        return(NULL)
    }

    data <- tryCatch(
        .build_association_plot_data(
            fit = x$fit,
            stan_data = x$stan_data,
            config = x$config,
            dataLong = x$dataLong,
            seed = seed
        ),
        error = function(e) NULL
    )
    if (!is.null(data)) {
        x$config$association_plot_data <- data
        x$cache_set("association_plot_data", data)
    }
    data
}

.association_support_range <- function(x, term_key, term = term_key, marker = NULL) {
    data <- .get_association_plot_data(x)
    support <- data$support
    if (is.null(support) || !nrow(support)) {
        return(NULL)
    }

    marker_key <- if (is.null(marker) || (length(marker) == 1L && is.na(marker))) "all" else as.character(marker[[1]])
    rows <- support[support$term == term & support$marker %in% c(marker_key, "all"), , drop = FALSE]
    if (!nrow(rows) && !identical(term, term_key)) {
        rows <- support[support$term == term_key & support$marker %in% c(marker_key, "all"), , drop = FALSE]
    }
    if (!nrow(rows)) {
        return(NULL)
    }
    c(min(rows$lower, na.rm = TRUE), max(rows$upper, na.rm = TRUE))
}

.build_association_plot_data <- function(fit, stan_data, config, dataLong, seed = 1) {
    assoc_flags <- config$assoc %||% list()
    active_terms <- names(assoc_flags)[vapply(assoc_flags, function(flag) isTRUE(as.logical(flag)), logical(1))]
    if (!length(active_terms)) {
        return(NULL)
    }

    data <- list(
        coeff_draws = list(),
        marker_weight_draws = list(),
        transform_coeff_draws = list(),
        transform_iota_draws = list(),
        transform_specs = config$transforms_spec %||% list(),
        support = data.frame(term = character(0), marker = character(0), lower = numeric(0), upper = numeric(0), source = character(0), stringsAsFactors = FALSE),
        term_map = data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE)
    )

    alpha_map <- c(
        cv_total = "alpha_cv_total",
        cv_mean = "alpha_cv_mean",
        cv_marker = "alpha_cv_marker",
        cs_total = "alpha_cs_total",
        cs_mean = "alpha_cs_mean",
        cs_marker = "alpha_cs_marker"
    )
    for (term_key in active_terms) {
        if (term_key %in% c("corr", "vcov")) {
            corr_vars <- grep(paste0("^alpha_", term_key, "_eff\\["), posterior::variables(.get_draws_obj(fit)), value = TRUE)
            if (!length(corr_vars)) {
                corr_vars <- grep(paste0("^alpha_", term_key, "\\["), posterior::variables(.get_draws_obj(fit)), value = TRUE)
            }
            if (length(corr_vars)) {
                data$coeff_draws[[term_key]] <- .get_draws_matrix(fit, variables = corr_vars, seed = seed)[, corr_vars, drop = FALSE]
                data$term_map <- rbind(
                    data$term_map,
                    data.frame(
                        term = sub("_eff\\[", "[", sub("^alpha_", "", corr_vars), perl = TRUE),
                        variable = corr_vars,
                        stringsAsFactors = FALSE
                    )
                )
            }
        } else {
            alpha_var <- alpha_map[[term_key]]
            if (!is.null(alpha_var)) {
                data$coeff_draws[[term_key]] <- .get_draws_matrix(fit, variables = alpha_var, seed = seed)[, 1, drop = FALSE]
                data$term_map <- rbind(data$term_map, data.frame(term = term_key, variable = alpha_var, stringsAsFactors = FALSE))
            }
        }

        tf_spec <- data$transform_specs[[term_key]] %||% list(type = "identity")
        tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")
        if (.is_ispline_transform_type(tf_type)) {
            map <- .assoc_channel_map(term_key)
            n_coeff <- as.integer(stan_data[[map$n_coeff]] %||% length(stan_data[[map$base_coeff]] %||% tf_spec$coeff %||% numeric(0)))
            if (n_coeff > 0L) {
                if (term_key %in% c("corr", "vcov")) {
                    n_components <- .assoc_transform_component_count(
                        term_key,
                        stan_data$Q_idm,
                        diagonal_only = identical(term_key, "vcov") && as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
                    )
                    component_terms <- .assoc_transform_component_labels(term_key, n_components)
                    for (m in seq_len(n_components)) {
                        eff_names <- paste0(map$eff_prefix, "[", m, ",", seq_len(n_coeff), "]")
                        dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                        if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                            data$transform_coeff_draws[[component_terms[[m]]]] <- as.matrix(dmat[, eff_names, drop = FALSE])
                        }
                    }
                } else {
                    eff_names <- paste0(map$eff_prefix, "[", seq_len(n_coeff), "]")
                    dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                    if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                        data$transform_coeff_draws[[term_key]] <- as.matrix(dmat[, eff_names, drop = FALSE])
                    }
                }
            }
        }

        if (identical(tf_type, "functional") && (isTRUE(tf_spec$estimate_iota_intercept) || isTRUE(tf_spec$estimate_iota_slope))) {
            map <- .assoc_channel_map(term_key)
            n_iota_intercept <- as.integer(tf_spec$n_iota_intercept %||% stan_data[[map$n_iota_intercept]] %||% 0L)
            n_iota_slope <- as.integer(tf_spec$n_iota_slope %||% stan_data[[map$n_iota_slope]] %||% 0L)

            get_cached_iota_draws <- function(prefix, n_iota, component_index = NULL) {
                if (is.null(prefix) || n_iota < 1L) {
                    return(NULL)
                }
                eff_names <- if (!is.null(component_index) && !is.na(component_index)) {
                    offsets <- (component_index - 1L) * n_iota + seq_len(n_iota)
                    paste0(prefix, "[", offsets, "]")
                } else {
                    paste0(prefix, "[", seq_len(n_iota), "]")
                }
                dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                    return(as.matrix(dmat[, eff_names, drop = FALSE]))
                }
                if (n_iota == 1L && is.null(component_index)) {
                    dmat_single <- tryCatch(.get_draws_matrix(fit, variables = prefix, seed = seed), error = function(e) NULL)
                    if (!is.null(dmat_single) && prefix %in% colnames(dmat_single)) {
                        return(matrix(as.numeric(dmat_single[, prefix, drop = TRUE]), ncol = 1L))
                    }
                }
                NULL
            }

            if (term_key %in% c("corr", "vcov")) {
                n_components <- .assoc_transform_component_count(
                    term_key,
                    stan_data$Q_idm,
                    diagonal_only = identical(term_key, "vcov") && as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
                )
                component_terms <- .assoc_transform_component_labels(term_key, n_components)
                for (m in seq_len(n_components)) {
                    data$transform_iota_draws[[component_terms[[m]]]] <- list(
                        intercept = get_cached_iota_draws(map$iota_intercept, n_iota_intercept, component_index = m),
                        slope = get_cached_iota_draws(map$iota_slope, n_iota_slope, component_index = m)
                    )
                }
            } else {
                data$transform_iota_draws[[term_key]] <- list(
                    intercept = get_cached_iota_draws(map$iota_intercept, n_iota_intercept),
                    slope = get_cached_iota_draws(map$iota_slope, n_iota_slope)
                )
            }
        }
    }

    if (any(active_terms %in% c("cv_total", "cs_total", "cv_marker", "cs_marker"))) {
        n_markers <- as.integer(stan_data$D %||% length(stan_data$marker_levels %||% numeric(0)))
        if (n_markers > 0L) {
            weight_terms <- intersect(active_terms, .weighted_assoc_term_keys())
            if (isTRUE(as.integer(stan_data$shared_marker_weights %||% 1L) == 1L) && length(weight_terms) > 1L) {
                weight_terms <- weight_terms[1L]
            }
            for (term_key in weight_terms) {
                eff_names <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(n_markers), "]")
                dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                if ((is.null(dmat) || !all(eff_names %in% colnames(dmat))) && isTRUE(as.integer(stan_data$shared_marker_weights %||% 1L) == 1L)) {
                    eff_names <- paste0("marker_weights_eff[", seq_len(n_markers), "]")
                    dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                }
                if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                    data$marker_weight_draws[[term_key]] <- as.matrix(dmat[, eff_names, drop = FALSE])
                }
            }
        }
    }

    data$support <- .model_implied_support(fit, stan_data, config, dataLong, seed = seed)
    data
}

.model_implied_support <- function(fit, stan_data, config, dataLong, seed = 1) {
    terms <- names(config$assoc %||% list())[vapply(config$assoc %||% list(), function(flag) isTRUE(as.logical(flag)), logical(1))]
    if (!length(terms) || is.null(stan_data$X_obs) || is.null(stan_data$id) || is.null(stan_data$marker)) {
        return(data.frame(term = character(0), marker = character(0), lower = numeric(0), upper = numeric(0), source = character(0), stringsAsFactors = FALSE))
    }

    var_names <- posterior::variables(.get_draws_obj(fit))
    mean_of_vars <- function(names_vec, default = 0) {
        if (!length(names_vec) || !all(names_vec %in% var_names)) {
            return(rep(default, length(names_vec)))
        }
        dmat <- .get_draws_matrix(fit, variables = names_vec, seed = seed)
        colMeans(dmat[, names_vec, drop = FALSE])
    }

    P <- as.integer(stan_data$P %||% ncol(stan_data$X_obs))
    beta_vars <- paste0("beta_eff_in_likelihood[", seq_len(P), "]")
    if (!all(beta_vars %in% var_names)) {
        beta_vars <- paste0("beta_scaled[", seq_len(P), "]")
    }
    if (!all(beta_vars %in% var_names)) {
        beta_vars <- paste0("beta[", seq_len(P), "]")
    }
    beta_mean <- mean_of_vars(beta_vars, default = 0)

    n_id <- as.integer(stan_data$n_id %||% max(stan_data$id))
    R_id <- as.integer(stan_data$R_id %||% ncol(stan_data$Z_id_obs) %||% 0L)
    R_mk <- as.integer(stan_data$R_mk %||% ncol(stan_data$Z_mk_obs) %||% 0L)
    Q_idm <- as.integer(stan_data$Q_idm %||% ncol(stan_data$Z_idm_obs) %||% 0L)
    D <- as.integer(stan_data$D %||% max(stan_data$marker))

    u_mean <- matrix(0, nrow = n_id, ncol = max(1L, R_id))
    if (R_id > 0L) {
        u_vars <- as.vector(outer(seq_len(n_id), seq_len(R_id), function(i, r) paste0("u_id[", i, ",", r, "]")))
        u_mean[] <- mean_of_vars(u_vars, default = 0)
    }
    v_mean <- matrix(0, nrow = max(1L, D), ncol = max(1L, R_mk))
    if (R_mk > 0L) {
        v_vars <- as.vector(outer(seq_len(D), seq_len(R_mk), function(d, r) paste0("v_marker[", d, ",", r, "]")))
        v_mean[] <- mean_of_vars(v_vars, default = 0)
    }
    w_mean <- array(0, dim = c(max(1L, n_id), max(1L, D), max(1L, Q_idm)))
    if (Q_idm > 0L) {
        w_vars <- character(0)
        for (i in seq_len(n_id)) {
            for (d in seq_len(D)) {
                for (q in seq_len(Q_idm)) {
                    w_vars <- c(w_vars, paste0("w_idm[", i, ",", d, ",", q, "]"))
                }
            }
        }
        vals <- mean_of_vars(w_vars, default = 0)
        dim(vals) <- c(n_id, D, Q_idm)
        w_mean <- vals
    }

    X_obs <- as.matrix(stan_data$X_obs)
    Z_id_obs <- as.matrix(stan_data$Z_id_obs %||% matrix(0, nrow(X_obs), R_id))
    Z_mk_obs <- as.matrix(stan_data$Z_mk_obs %||% matrix(0, nrow(X_obs), R_mk))
    Z_idm_obs <- as.matrix(stan_data$Z_idm_obs %||% matrix(0, nrow(X_obs), Q_idm))
    id_idx <- as.integer(stan_data$id)
    marker_idx <- as.integer(stan_data$marker)

    cv_mean <- as.numeric(X_obs %*% beta_mean)
    if (R_id > 0L) {
        cv_mean <- cv_mean + rowSums(Z_id_obs * u_mean[id_idx, seq_len(R_id), drop = FALSE])
    }

    cv_marker <- rep(0, nrow(X_obs))
    if (R_mk > 0L) {
        cv_marker <- cv_marker + rowSums(Z_mk_obs * v_mean[marker_idx, seq_len(R_mk), drop = FALSE])
    }
    if (Q_idm > 0L) {
        w_rows <- do.call(rbind, lapply(seq_len(nrow(X_obs)), function(n) w_mean[id_idx[n], marker_idx[n], seq_len(Q_idm)]))
        cv_marker <- cv_marker + rowSums(Z_idm_obs * w_rows)
    }
    cv_total <- cv_mean + cv_marker

    marker_var <- if ("marker" %in% names(dataLong)) "marker" else names(dataLong)[match(TRUE, grepl("marker", names(dataLong), ignore.case = TRUE))]
    if (is.na(marker_var) || is.null(marker_var)) marker_var <- "marker"
    marker_labels <- as.character(stan_data$marker_levels %||% sort(unique(as.character(dataLong[[marker_var]]))))
    if (!length(marker_labels)) {
        marker_labels <- as.character(sort(unique(marker_idx)))
    }
    row_markers <- if (!is.null(stan_data$marker_levels)) as.character(stan_data$marker_levels[marker_idx]) else as.character(marker_idx)
    time_var <- if ("time" %in% names(dataLong)) "time" else names(dataLong)[match(TRUE, grepl("time", names(dataLong), ignore.case = TRUE))]
    if (is.na(time_var) || is.null(time_var) || !(time_var %in% names(dataLong))) {
        time_vals <- seq_len(nrow(dataLong))
    } else {
        time_vals <- dataLong[[time_var]]
    }

    make_support_rows <- function(term_key, values, source = "model_implied") {
        rows <- list()
        idx <- 1L
        all_q <- stats::quantile(values, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
        if (all(is.finite(all_q)) && diff(all_q) > 0) {
            rows[[idx]] <- data.frame(term = term_key, marker = "all", lower = all_q[1], upper = all_q[2], source = source, stringsAsFactors = FALSE)
            idx <- idx + 1L
        }
        for (mk in unique(row_markers)) {
            mk_vals <- values[row_markers == mk]
            mk_q <- stats::quantile(mk_vals, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
            if (all(is.finite(mk_q)) && diff(mk_q) > 0) {
                rows[[idx]] <- data.frame(term = term_key, marker = as.character(mk), lower = mk_q[1], upper = mk_q[2], source = source, stringsAsFactors = FALSE)
                idx <- idx + 1L
            }
        }
        do.call(rbind, rows)
    }

    slope_by_group <- function(values) {
        grp <- interaction(id_idx, row_markers, drop = TRUE)
        out <- rep(NA_real_, length(values))
        for (g in unique(grp)) {
            ii <- which(grp == g)
            if (length(ii) < 2L) next
            ord <- order(time_vals[ii])
            dt <- diff(as.numeric(time_vals[ii][ord]))
            dv <- diff(values[ii][ord])
            sl <- dv / dt
            out[ii[ord][-1]] <- sl
        }
        out[is.finite(out)]
    }

    out <- list()
    if ("cv_mean" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cv_mean", cv_mean)
    if ("cv_marker" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cv_marker", cv_marker)
    if ("cv_total" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cv_total", cv_total)
    if ("cs_mean" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cs_mean", slope_by_group(cv_mean))
    if ("cs_marker" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cs_marker", slope_by_group(cv_marker))
    if ("cs_total" %in% terms) out[[length(out) + 1L]] <- make_support_rows("cs_total", slope_by_group(cv_total))
    if ("corr" %in% terms) {
        n_corr <- .assoc_transform_component_count("corr", as.integer(stan_data$Q_idm), diagonal_only = FALSE)
        corr_terms <- .assoc_transform_component_labels("corr", n_corr)
        out[[length(out) + 1L]] <- do.call(rbind, lapply(corr_terms, function(term_label) {
            data.frame(term = term_label, marker = "all", lower = -1, upper = 1, source = "theoretical", stringsAsFactors = FALSE)
        }))
    }
    if ("vcov" %in% terms && isTRUE((stan_data$Q_idm %||% 0L) > 0L)) {
        vcov_map <- .assoc_cov_feature_map(
            q_idm = as.integer(stan_data$Q_idm),
            diagonal_only = as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L,
            include_diag = TRUE
        )
        alpha_mean <- mean_of_vars(paste0("alpha_L[", seq_len(nrow(vcov_map)), "]"), default = 0)
        lambda_mean <- mean_of_vars(paste0("lambda_L[", seq_len(nrow(vcov_map)), "]"), default = 0)
        z_l_names <- as.vector(outer(seq_len(n_id), seq_len(nrow(vcov_map)), function(i, m) paste0("z_L[", i, ",", m, "]")))
        z_l_mean <- if (length(z_l_names) > 0L) {
            matrix(
                mean_of_vars(z_l_names, default = 0),
                nrow = n_id,
                ncol = nrow(vcov_map),
                byrow = TRUE
            )
        } else {
            matrix(0, nrow = n_id, ncol = 0L)
        }
        k_cov <- as.integer(stan_data$K_cov %||% 0L)
        beta_mean <- if (k_cov > 0L && nrow(vcov_map) > 0L) {
            beta_flat <- mean_of_vars(as.vector(outer(seq_len(nrow(vcov_map)), seq_len(k_cov), function(m, k) paste0("beta_L[", m, ",", k, "]"))), default = 0)
            matrix(beta_flat, nrow = nrow(vcov_map), byrow = TRUE)
        } else {
            matrix(0, nrow = nrow(vcov_map), ncol = 0L)
        }
        xcov <- as.matrix(stan_data$Xcov %||% matrix(0, nrow = n_id, ncol = k_cov))
        marker_id_row_scale <- rep(1, as.integer(stan_data$Q_idm %||% 0L))
        idx_time_idm <- as.integer(stan_data$idx_time_idm %||% integer(0))
        idx_time_idm <- idx_time_idm[is.finite(idx_time_idm) & idx_time_idm >= 1L & idx_time_idm <= length(marker_id_row_scale)]
        if (length(idx_time_idm) > 0L) {
            marker_id_row_scale[idx_time_idm] <- as.numeric(stan_data$tmax_internal %||% 1)
        }
        li_terms <- lapply(seq_len(n_id), function(i) {
            if (nrow(vcov_map) == 0L) {
                return(numeric(0))
            }
            lp_vec <- vapply(seq_len(nrow(vcov_map)), function(m) {
                alpha_mean[m] +
                    if (k_cov > 0L) sum(beta_mean[m, ] * xcov[i, ]) else 0 +
                    lambda_mean[m] * z_l_mean[min(i, nrow(z_l_mean)), m]
            }, numeric(1))
            li <- .cov_lp_to_chol(
                lp_vec = lp_vec,
                q_idm = as.integer(stan_data$Q_idm),
                idx_row = vcov_map[, 1],
                idx_col = vcov_map[, 2],
                diag_link = stan_data$vcov_diag_link
            )
            li <- sweep(li, 1L, marker_id_row_scale, `*`)
            .assoc_vcov_features_from_chol(
                li,
                diagonal_only = as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
            )
        })
        if (length(li_terms) > 0L && length(li_terms[[1]]) > 0L) {
            li_mat <- do.call(rbind, li_terms)
            term_labels <- .assoc_transform_component_labels("vcov", ncol(li_mat))
            for (j in seq_len(ncol(li_mat))) {
                vals_j <- li_mat[, j]
                vals_j <- vals_j[is.finite(vals_j)]
                if (length(vals_j) > 1L) {
                    out[[length(out) + 1L]] <- make_support_rows(term_labels[[j]], vals_j)
                }
            }
        }
    }

    out <- out[!vapply(out, is.null, logical(1))]
    if (!length(out)) {
        return(data.frame(term = character(0), marker = character(0), lower = numeric(0), upper = numeric(0), source = character(0), stringsAsFactors = FALSE))
    }
    do.call(rbind, out)
}

.joinme_make_assoc_transform <- function(spec, term_name) {
    if (is.null(spec) || is.null(spec$type) || identical(spec$type, "identity")) {
        return(function(x) x)
    }

    tf_type <- .canonicalise_transform_type(spec$type)

    if (identical(tf_type, "functional")) {
        bc <- parse_transform_expr(spec$expr, iota_nodes = .transform_iota_nodes(spec))
        return(function(x, iota_intercept = numeric(0), iota_slope = numeric(0)) {
            eval_bytecode_vector(
                x = x,
                bytecode = bc$bytecode %||% bc$opcodes,
                const_data = bc$const_data %||% numeric(0),
                iota_intercepts = iota_intercept,
                iota_slopes = iota_slope,
                op_iota_intercept_idx = bc$op_iota_intercept_idx %||% integer(length(bc$bytecode %||% bc$opcodes %||% integer(0))),
                op_iota_slope_idx = bc$op_iota_slope_idx %||% integer(length(bc$bytecode %||% bc$opcodes %||% integer(0)))
            )
        })
    }

    if (.is_ispline_transform_type(tf_type)) {
        if (tf_type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised") && is.null(spec$y)) {
            spec <- .make_stan_penalised_ispline_transform(spec)
        } else if (tf_type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised")) {
            spec <- .make_penalised_ispline_transform(spec)
        }
        if (!requireNamespace("splines2", quietly = TRUE)) {
            cli::cli_abort(c(
                x = "Transform type {.val {tf_type}} for {.val {term_name}} requires {.pkg splines2}.",
                i = "Install {.pkg splines2} or use identity/functional/pwlin transforms."
            ))
        }
        knots <- as.numeric(spec$knots)
        if (.transform_uses_expit_input(spec)) {
            knots <- .validate_expit_domain_values(knots, "knots")
        }
        coeff <- as.numeric(spec$coeff)
        degree <- as.integer(spec$degree %||% 3L)
        return(function(x) {
            x <- .transform_input_for_spec(x, spec)
            boundary <- c(knots[1], knots[length(knots)])
            x <- pmin(pmax(x, boundary[1]), boundary[2])
            internal_knots <- if (length(knots) > 2) knots[2:(length(knots) - 1)] else numeric(0)
            basis <- splines2::iSpline(x,
                knots = internal_knots,
                degree = degree,
                intercept = TRUE,
                Boundary.knots = boundary
            )
            coef_len <- ncol(basis)
            coeff_use <- if (length(coeff) < coef_len) c(coeff, rep(0, coef_len - length(coeff))) else coeff[seq_len(coef_len)]
            vals <- as.numeric(basis %*% coeff_use)
            if (.resolve_monotone_direction(spec$direction %||% spec$spline_direction %||% 1L, default = 1L) < 0L) {
                vals <- sum(coeff_use) - vals
            }
            vals
        })
    }

    if (identical(tf_type, "pwlin")) {
        xk <- as.numeric(spec$x)
        yk <- as.numeric(spec$y)
        return(function(x) stats::approx(x = xk, y = yk, xout = as.numeric(x), method = "linear", rule = 2)$y)
    }

    cli::cli_abort(c(
        x = "Unsupported transform type for {.val {term_name}}: {.val {tf_type}}.",
        i = "Use one of identity, functional, ispline, ispline_penalised, ispline_expit, ispline_expit_penalised, or pwlin."
    ))
}

.running_mean_df <- function(arr) {
    dims <- dim(arr)
    n_iter <- dims[1]
    n_chain <- dims[2]
    n_var <- dims[3]
    var_names <- dimnames(arr)[[3]]
    out <- vector("list", n_var * n_chain)
    idx <- 1
    for (v in seq_len(n_var)) {
        for (ch in seq_len(n_chain)) {
            vec <- arr[, ch, v]
            run_mean <- cumsum(vec) / seq_along(vec)
            out[[idx]] <- data.frame(
                iteration = seq_len(n_iter),
                chain = ch,
                variable = var_names[[v]],
                value = run_mean,
                stringsAsFactors = FALSE
            )
            idx <- idx + 1
        }
    }
    do.call(rbind, out)
}

.running_quantile_df <- function(arr, probs = c(0.1, 0.5, 0.9)) {
    dims <- dim(arr)
    n_iter <- dims[1]
    n_chain <- dims[2]
    n_var <- dims[3]
    var_names <- dimnames(arr)[[3]]
    stat_names <- vapply(probs, .quantile_name_from_prob, character(1))

    out <- vector("list", n_var * n_chain)
    idx <- 1
    for (v in seq_len(n_var)) {
        for (ch in seq_len(n_chain)) {
            vec <- arr[, ch, v]
            qmat <- vapply(seq_len(n_iter), function(i) {
                stats::quantile(vec[seq_len(i)], probs = probs, names = FALSE)
            }, numeric(length(probs)))
            df <- data.frame(
                iteration = rep(seq_len(n_iter), each = length(probs)),
                chain = ch,
                variable = var_names[[v]],
                stat = rep(stat_names, times = n_iter),
                value = as.vector(qmat),
                stringsAsFactors = FALSE
            )
            out[[idx]] <- df
            idx <- idx + 1
        }
    }
    do.call(rbind, out)
}

# ============================================================================
# Utility: NULL or default operator
# ============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x
