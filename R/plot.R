#' Enhanced Plot Method for Dynamic Prediction from Joint Models
#'
#' @importFrom stats na.omit
#' @importFrom dplyr %>% all_of
#'
#' @description
#' Produces publication-ready ggplot2 visualizations of dynamic predictions from
#' a joint model with flexible customization of longitudinal trajectories and
#' conditional survival curves.
#'
#' @param x An object of class `JoinMeDynPred` returned by the [predict] method.
#' @param which Character vector indicating what to plot: "cumhaz" (default),
#'   "longitudinal", "survival", or "both".
#' @param subject Integer/Character vector. Which subject(s) to plot. If NULL,
#'   plots all subjects in separate panels or returns a list.
#' @param trajectory_type Character. Which trajectory to show:
#'   "id_marker" (default), "marker_pop" (marker-level average),
#'   "overall_pop" (overall mean), or "all_types" for all three.
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
#' @param combined Logical. If TRUE and both outcomes requested, create a combined plot.
#'   Otherwise return separate plots.
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
#' - The longitudinal scale is inferred from the prediction object.
#' - On \"predict\" scale: observed values are shown for times in interval
#'   `[0, T_start]` and predictions are drawn for `[T_start, T_horiz]`.
#' - On \"linpred\" or \"epred\" scale: predictions are shown for `[0, T_horiz]`,
#'   with fitted values providing the pre-T_start segment when available.
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
#' If multiple subjects: a list of plots or a combined faceted plot.
#' If both outcomes: list with "longitudinal", "survival", and "cumhaz" elements,
#'   or a single combined plot if `combined = TRUE`.
#'
#' @export
plot.JoinMeDynPred <- function(
    x,
    which = c("cumhaz", "longitudinal", "survival", "both"),
    subject = NULL,
    trajectory_type = c("id_marker", "marker_pop", "overall_pop", "all_types"),
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
    is_competing <- FALSE
    if (!is.null(x$data$event)) {
        ev_cols <- names(x$data$event)
        is_competing <- any(c("event_type", "stage_from", "stage_to") %in% ev_cols)
    }
    if (missing(which) && isTRUE(is_competing)) {
        which <- "cumhaz"
    }
    which <- match.arg(which)
    trajectory_type <- match.arg(trajectory_type)
    smooth_method <- match.arg(smooth_method)
    facet_by <- match.arg(facet_by)
    ci_type <- match.arg(ci_type)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    ci_levels <- .validate_ci_levels_plot(ci_levels, x$metadata$ci_levels %||% NULL)

    # Check for required data
    if (is.null(x$quantiles)) {
        stop("This version of plot() requires enhanced prediction output with quantiles.\n",
             "Please update your predict() output or re-run prediction.")
    }

    # Identify subjects
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
        if (which %in% c("both", "longitudinal")) {
            for (traj_src in trajectory_sources) {
                p_long <- .plot_longitudinal_single(
                    x, id, traj_src,
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
        if (which %in% c("both", "survival")) {
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
        if (which %in% c("both", "cumhaz")) {
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
        } else if (isTRUE(combined) && !is.null(res$longitudinal) && !is.null(res$survival)) {
            # Combine longitudinal and survival (cumhaz returned separately if present)
            combined_plot <- .combine_long_surv_plots(res$longitudinal, res$survival)
            if (!is.null(res$cumhaz)) {
                return(list(combined = combined_plot, cumhaz = res$cumhaz))
            }
            return(combined_plot)
        } else {
            return(res)
        }
    } else {
        # Multiple subjects
        if (isTRUE(combined) && which == "both") {
            # Create a grid with all subjects
            all_p_long <- lapply(plots_list, function(x) x$longitudinal)
            all_p_surv <- lapply(plots_list, function(x) x$survival)
            all_p_cumhaz <- lapply(plots_list, function(x) x$cumhaz)
            return(list(
                longitudinal = all_p_long,
                survival = all_p_surv,
                cumhaz = all_p_cumhaz
            ))
        } else if (which == "longitudinal") {
            return(lapply(plots_list, function(x) x$longitudinal))
        } else if (which == "survival") {
            return(lapply(plots_list, function(x) x$survival))
        } else if (which == "cumhaz") {
            return(lapply(plots_list, function(x) x$cumhaz))
        }
        return(plots_list)
    }
}

# ============================================================================
# Helper: Longitudinal Plot for Single Subject
# ============================================================================
.plot_longitudinal_single <- function(
    x, id, trajectory_type = "id_marker",
    smooth_trajectory, smooth_method, smooth_span,
    ci_levels, ci_type, observed_first,
    facet_by, facet_scales,
    show_data, show_observed_line,
    observed_style, prediction_style,
    theme_fn, palette_marker
) {

    # Extract quantile data based on trajectory type
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
    
    if (is.null(quant_df)) quant_df <- x$predictions$longitudinal
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
    
    if (nrow(quant_df) == 0) {
        cli::cli_warn(c(
            x = "No data found for subject {id} ({trajectory_type}).",
            i = "Check the prediction output and subject identifiers."
        ))
        return(NULL)
    }
    scale_long <- x$metadata$scale %||% "epred"
    if (!scale_long %in% c("epred", "linpred", "predict")) scale_long <- "epred"

    # Get observed data
    data_long <- x$data$longitudinal
    data_long_id <- data_long[data_long$id == id, , drop = FALSE]

    # Identify variable names from metadata (stored during prediction)
    # Fall back to x$call if not available
    id_var <- x$metadata$id_var %||% eval(x$call$id_var) %||% "id"
    time_var <- x$metadata$time_var %||% eval(x$call$time_var) %||% "time"
    marker_var <- x$metadata$marker_var %||% eval(x$call$marker_var) %||% "marker"
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
            stop("Cannot determine response variable name from prediction object. ",
                 "Please ensure the prediction object was created with an updated version of predict().")
        }
    }

    # Conditioning time
    t_cond <- x$metadata$conditioning_time

    # Initialize plot data
    df_pred <- quant_df |>
        dplyr::mutate(
            type = "prediction",
            marker = as.character(.data$marker)
        )

    # Prepare observed or fitted data for the left-of-Tstart region
    df_obs <- NULL
    use_fitted <- scale_long %in% c("linpred", "epred")

    if (isTRUE(show_data)) {
        if (use_fitted) {
            fit_df <- x$quantiles$longitudinal_fitted
            if (!is.null(fit_df) && "scale" %in% names(fit_df)) {
                fit_df <- fit_df[fit_df$id == id & fit_df$scale == scale_long, , drop = FALSE]
                if (nrow(fit_df) > 0) {
                    median_col <- .quantile_name_from_prob(0.5)
                    if (!(median_col %in% names(fit_df))) median_col <- "mean"
                    df_obs <- fit_df |>
                        dplyr::select(time, marker, all_of(median_col)) |>
                        dplyr::rename(value = all_of(median_col)) |>
                        dplyr::mutate(
                            type = ifelse(.data$time <= t_cond, "fitted", "fitted_future"),
                            marker = as.character(.data$marker)
                        )
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
                df_obs <- data_long_id |>
                    dplyr::select(all_of(c(time_var, marker_var, resp_var))) |>
                    dplyr::rename(time = all_of(time_var), marker = all_of(marker_var)) |>
                    dplyr::mutate(
                        type = ifelse(.data$time <= t_cond, "observed", "observed_future"),
                        marker = as.character(.data$marker)
                    ) |>
                    dplyr::rename(value = all_of(resp_var))
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
        df_pred_future <- df_pred[df_pred$time <= t_horiz, , drop = FALSE]
        df_pred_hist <- NULL
        fit_quant <- x$quantiles$longitudinal_fitted
        if (!is.null(fit_quant)) {
            fit_quant <- fit_quant[fit_quant$id == id & fit_quant$scale == scale_long, , drop = FALSE]
            if (nrow(fit_quant) > 0) {
                df_pred_hist <- fit_quant
            }
        }
        df_pred <- if (!is.null(df_pred_hist)) {
            dplyr::bind_rows(df_pred_hist, df_pred_future)
        } else {
            df_pred_future
        }
    }

    # Start plot with data and base aesthetics
    if (is.null(df_pred) || nrow(df_pred) == 0) {
        stop("No valid data to plot for subject ", id)
    }
    
    # Initialize plot - don't set y in base aes since different layers use different y columns
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
                df_ribbon <- df_pred |>
                    dplyr::select(time, marker, all_of(p_names)) |>
                    dplyr::rename(ymin = all_of(p_names[1]), ymax = all_of(p_names[2]))

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
                    df_line <- df_pred |> dplyr::select(time, marker, all_of(qcol)) |>
                        dplyr::rename(value = all_of(qcol))
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
    if (nrow(quant_df) == 0) {
        cli::cli_warn(c(
            x = "No survival data found for subject {id}.",
            i = "Ensure survival predictions were requested."
        ))
        return(NULL)
    }

    # Initialize plot with data and base aesthetics (use q50 as y for survival)
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
                # Use aes_string for safe dynamic column references in loops
                p <- p + ggplot2::geom_ribbon(
                    ggplot2::aes_string(
                        x = "time",
                        ymin = p_names[1],
                        ymax = p_names[2]
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
                        ggplot2::aes_string(x = "time", y = qcol),
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
                    ggplot2::aes_string(
                        x = "time",
                        ymin = p_names[1],
                        ymax = p_names[2]
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
                        ggplot2::aes_string(x = "time", y = qcol),
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
    lbl <- sub("\\.?0+$", "", lbl)
    paste0("q", lbl)
}

.quantile_names_from_ci <- function(level) {
    .quantile_name_from_prob(c((1 - level) / 2, (1 + level) / 2))
}

.latex_label_long <- function(resp_var, scale_label) {
    base_label <- paste0(resp_var, " (", scale_label, ")")
    if (!is.na(resp_var) && requireNamespace("latex2exp", quietly = TRUE)) {
        latex2exp::TeX(paste0("$", resp_var, "$\\,(\\mathrm{", scale_label, "})"))
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
# Utility: NULL or default operator
# ============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x
