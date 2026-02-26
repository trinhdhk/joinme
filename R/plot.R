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
#'   "longitudinal", or "survival". Multiple values are allowed.
#' @param subject Integer/Character vector. Which subject(s) to plot. If NULL,
#'   plots all subjects in separate panels or returns a list.
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
    which = c("cumhaz", "longitudinal", "survival"),
    subject = NULL,
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
    is_competing <- FALSE
    if (!is.null(x$data$event)) {
        ev_cols <- names(x$data$event)
        is_competing <- any(c("event_type", "stage_from", "stage_to") %in% ev_cols)
    }
    if (missing(which)) {
        which <- if (isTRUE(is_competing)) "cumhaz" else c("longitudinal", "survival")
    }
    which <- .validate_plot_which(which)
    trajectory_type <- match.arg(trajectory_type)
    smooth_method <- match.arg(smooth_method)
    facet_by <- match.arg(facet_by)
    ci_type <- match.arg(ci_type)

    ci_levels <- .validate_ci_levels_plot(ci_levels, x$metadata$ci_levels %||% NULL)

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
        if ("longitudinal" %in% which) {
            for (traj_src in trajectory_sources) {
                p_long <- .plot_longitudinal_single(
                    x, id, traj_src,
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
        if ("survival" %in% which) {
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
        if ("cumhaz" %in% which) {
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
                combined_plot <- .combine_plot_grid(res, ncol = length(which))
                if (inherits(combined_plot, "gg")) return(combined_plot)
            }
            return(res)
        }

        # Multiple subjects
        if (isTRUE(combined)) {
            # Combine requested outcomes within each subject.
            # Output contract for multiple subjects:
            # - preferred: list(id -> combined plot object)
            # - fallback: list(id -> uncombined per-subject plot list)
            out <- lapply(plots_list, function(sub_plots) {
                design <- switch(
                    length(which),
                    "1" = 'A',
                    "2" = "AB",
                    "3" = "AB\nAC"
                )
                .combine_plot_grid(
                    sub_plots,
                    ncol = min(2, length(which)),
                    nrow = min(1, length(which) - 1),
                    design = design,
                    fallback = "input"
                )
            })
            return(out)
        }
        if (length(which) == 1 && which == "longitudinal") {
            return(lapply(plots_list, function(x) x$longitudinal))
        }
        if (length(which) == 1 && which == "survival") {
            return(lapply(plots_list, function(x) x$survival))
        }
        if (length(which) == 1 && which == "cumhaz") {
            return(lapply(plots_list, function(x) x$cumhaz))
        }
        plots_list
}

# ============================================================================
# Helper: Longitudinal Plot for Single Subject
# ============================================================================
.plot_longitudinal_single <- function(
    x, id, trajectory_type = "id_marker",
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

    # Initialize plot data
    df_pred <- quant_df |>
        dplyr::mutate(
            type = "prediction",
            segment = "predict",
            marker = as.character(.data$marker)
        )

    # Prepare observed or fitted data for the left-of-time_start region
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
        df_pred_future <- df_pred[df_pred$time >= t_cond & df_pred$time <= t_horiz, , drop = FALSE]
        df_pred_future$segment <- "predict"
        df_pred_hist <- NULL
        fit_quant <- x$quantiles$longitudinal_fitted
        if (!is.null(fit_quant)) {
            fit_quant <- fit_quant[fit_quant$id == id & fit_quant$scale == scale_long, , drop = FALSE]
            if (nrow(fit_quant) > 0) {
                df_pred_hist <- fit_quant
                df_pred_hist$segment <- "fitted"
                df_pred_hist <- df_pred_hist[df_pred_hist$time <= t_cond, , drop = FALSE]
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
.validate_plot_which <- function(which) {
    allowed <- c("cumhaz", "longitudinal", "survival")
    which <- unique(which)
    bad <- setdiff(which, allowed)
    if (length(bad) > 0) {
        cli::cli_abort(c(
            x = "Unknown plot type(s): {paste(bad, collapse = ', ')}.",
            i = "Use one or more of: {paste(allowed, collapse = ', ')}."
        ))
    }
    which
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

#' Diagnostic plots for JoinMe model fits
#'
#' @param x A fitted object of class `JoinMeFit`.
#' @param type Diagnostic plot type.
#' @param pars Optional character vector of parameter names to include.
#' @param regex_pars Optional regular expression for parameter selection.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param max_vars Maximum number of parameters for running diagnostics.
#' @param quantile_probs Numeric vector of quantile probabilities for running quantile plots.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
plot.JoinMeFit <- function(x, type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean", "running_quantile"),
                           pars = NULL, regex_pars = NULL, draws = NULL, seed = 1, max_vars = 4,
                           quantile_probs = c(0.1, 0.5, 0.9), ...) {
    if (!inherits(x, "JoinMeFit")) {
        cli::cli_abort("{.arg x} must be a JoinMeFit object.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }
    type <- match.arg(type)

    if (type %in% c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd")) {
        df <- switch(type,
            rhat = stan_rhat(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed),
            ess_bulk = stan_ess(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "bulk"),
            ess_tail = stan_ess(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "tail"),
            mcse_mean = stan_mcse(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "mean"),
            mcse_sd = stan_mcse(x, pars = pars, regex_pars = regex_pars, draws = draws, seed = seed, type = "sd")
        )
        if (nrow(df) == 0) {
            cli::cli_abort("No parameters found for diagnostic plotting.")
        }
        metric_name <- setdiff(names(df), "variable")
        names(df)[names(df) == metric_name] <- "metric"
        p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(variable, metric), y = metric)) +
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
        p <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = value, color = factor(chain))) +
            ggplot2::geom_line(alpha = 0.7) +
            ggplot2::facet_wrap(~ variable, scales = "free_y") +
            ggplot2::labs(x = "Iteration", y = "Running mean", color = "Chain", title = "Running mean by chain") +
            ggplot2::theme_minimal()
        return(p)
    }

    df <- .running_quantile_df(arr, probs = quantile_probs)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = value, color = factor(chain))) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::facet_grid(variable ~ stat, scales = "free_y") +
        ggplot2::labs(x = "Iteration", y = "Running quantile", color = "Chain", title = "Running quantiles by chain") +
        ggplot2::theme_minimal()
    p
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
