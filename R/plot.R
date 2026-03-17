#' Enhanced Plot Method for Dynamic Prediction from Joint Models
#'
#' @importFrom stats na.omit
#' @importFrom dplyr %>% all_of
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
            # Output contract for multiple subjects:
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
                if (!is.null(marker) && "marker" %in% names(fit_df)) {
                    fit_df <- fit_df[as.character(fit_df$marker) %in% marker, , drop = FALSE]
                }
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
            dplyr::bind_rows(df_pred_hist, df_pred_future)
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
#'   `"survival"`, `"cumhaz"`, and `"association"`.
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
#'   $\alpha$; `association_metric = "transform"` plots the transform $f(x)$
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
                           conditioning = c("auto", "last", "origin"),
                           time_start = NULL,
                           times = NULL,
                           time_horizon = NULL,
                           pred_control = list(n_samples = 100),
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
    fitted_types <- c("longitudinal", "survival", "cumhaz", "association")

    if (!inherits(x, "JoinMeFit")) {
        cli::cli_abort("{.arg x} must be a JoinMeFit object.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }

    if ("which" %in% names(dots)) {
        type <- dots$which
    }
    type <- unique(as.character(type))
    bad_types <- setdiff(type, c(diagnostic_types, fitted_types))
    if (length(bad_types) > 0) {
        cli::cli_abort(c(
            x = "Unknown {.arg type}: {paste(bad_types, collapse = ', ')}.",
            i = "Use one or more of: {paste(c(diagnostic_types, fitted_types), collapse = ', ')}."
        ))
    }

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

    conditioning <- match.arg(conditioning)
    smooth_method <- match.arg(smooth_method)
    facet_by <- match.arg(facet_by)
    ci_type <- match.arg(ci_type)
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
        conditioning = conditioning,
        time_start = time_start,
        times = times,
        time_horizon = time_horizon,
        pred_control = pred_control,
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
            i = "Supported names are association_term, association_grid, association_range, association_points, and association_metric."
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

    legacy_names <- names(opts)
    for (nm in intersect(legacy_names, names(dots))) {
        if (nm %in% names(association_options)) {
            next
        }
        opts[[nm]] <- dots[[nm]]
    }

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

.plot_joinmefit_fitted <- function(x, type, subject, marker, scale, conditioning, time_start, times, time_horizon,
                                   pred_control, seed, smooth_trajectory, smooth_method, smooth_span,
                                   ci_levels, ci_type, observed_first, facet_by, facet_scales, combined,
                                   show_data, show_observed_line, observed_style, prediction_style,
                                   theme_fn, palette_marker) {
    pred <- .build_joinmefit_plot_prediction(
        x = x,
    which = type,
        subject = subject,
        scale = scale,
        conditioning = conditioning,
        time_start = time_start,
        times = times,
        time_horizon = time_horizon,
        pred_control = pred_control,
        seed = seed,
        ci_levels = ci_levels
    )

    plot(pred,
            type = intersect(type, c("longitudinal", "survival", "cumhaz")),
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

.build_joinmefit_plot_prediction <- function(x, which, subject, scale, conditioning, time_start, times,
                                             time_horizon, pred_control, seed, ci_levels) {
    id_var <- .joinmefit_call_arg_chr(x$call, "id_var", "id")
    time_var <- .joinmefit_call_arg_chr(x$call, "time_var", "time")

    data_long <- x$dataLong
    data_event <- x$dataEvent
    if (!is.null(subject)) {
        keep_ids <- as.character(subject)
        data_long <- data_long[as.character(data_long[[id_var]]) %in% keep_ids, , drop = FALSE]
        data_event <- data_event[as.character(data_event[[id_var]]) %in% keep_ids, , drop = FALSE]
    }

    time_start_use <- .resolve_joinmefit_plot_time_start(
        data_long = data_long,
        id_var = id_var,
        time_var = time_var,
        which = which,
        conditioning = conditioning,
        time_start = time_start
    )

    process <- c(
        if ("longitudinal" %in% which) "longitudinal",
        if ("longitudinal" %in% which || any(c("survival", "cumhaz") %in% which)) "event"
    )
    process <- unique(process)

    scale_use <- if ("longitudinal" %in% which) .normalize_prediction_scales(scale %||% c("epred", "linpred", "predict")) else NULL
    control_use <- utils::modifyList(list(n_samples = 100L), pred_control %||% list())

    predict.JoinMeFit(
        object = x,
        newdataLong = data_long,
        newdataEvent = data_event,
        process = process,
        pred_type = "per_marker_id",
        scale = scale_use,
        times = times,
        time_start = time_start_use,
        time_horizon = time_horizon,
        ci_levels = ci_levels,
        control = control_use,
        seed = seed
    )
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
    assoc_terms <- .joinmefit_available_association_terms(x)
    if (!is.null(association_term)) {
        assoc_terms <- .joinmefit_expand_association_terms(association_term, assoc_terms)
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

.joinmefit_expand_association_terms <- function(requested_terms, available_terms) {
    requested_terms <- unique(as.character(requested_terms %||% character(0)))
    available_terms <- unique(as.character(available_terms %||% character(0)))
    if (!length(requested_terms) || !length(available_terms)) {
        return(character(0))
    }

    matched <- unlist(lapply(requested_terms, function(requested_term) {
        hits <- available_terms[
            available_terms == requested_term |
                startsWith(available_terms, paste0(requested_term, "["))
        ]
        hits
    }), use.names = FALSE)

    unique(matched)
}

.joinmefit_available_association_terms <- function(x) {
    payload <- .joinmefit_get_association_plot_payload(x)
    if (!is.null(payload$term_map) && nrow(payload$term_map) > 0L) {
        terms <- unique(as.character(payload$term_map$term))
        return(terms[!grepl("^weight:\\s", terms)])
    }
    assoc_draws <- extract.JoinMeFit(x, what = "assoc", keep_chains = FALSE)
    terms <- unique(as.character(assoc_draws$term_map$term))
    terms[!grepl("^weight:\\s", terms)]
}

.plot_joinmefit_association_single <- function(x, term, marker, association_grid, association_range, association_points, ci_levels,
                                               ci_type, association_metric, prediction_style, theme_fn, seed) {
    term_key <- if (grepl("^corr", term)) {
        "corr"
    } else if (grepl("^vcov", term)) {
        "vcov"
    } else {
        term
    }
    payload <- .joinmefit_get_association_plot_payload(x, seed = seed)
    coeff_draws <- .joinmefit_assoc_coeff_draws(x, term = term, payload = payload, seed = seed)
    marker_levels <- .association_plot_markers(x, term_key, marker)
    probs <- .quantile_probs_from_ci_plot(ci_levels)
    q_names <- .quantile_colnames(probs)
    plot_df <- do.call(rbind, lapply(marker_levels, function(marker_level) {
        x_grid <- association_grid %||% .default_joinmefit_association_grid(
            x,
            term_key,
            marker = marker_level,
            association_range = association_range,
            n = association_points
        )
        x_grid <- sort(unique(as.numeric(x_grid)))
        tf_mat <- .joinmefit_association_transform_matrix(x, term_key, term = term, x_grid = x_grid, n_draws = length(coeff_draws), seed = seed, payload = payload)
        weight_draws <- .joinmefit_association_marker_weight_draws(
            x = x,
            term_key = term_key,
            marker_level = marker_level,
            n_draws = length(coeff_draws),
            seed = seed,
            payload = payload
        )
        tf_mat <- tf_mat * as.numeric(weight_draws)

        if (identical(association_metric, "hazard") && term_key %in% c("corr", "vcov")) {
            tf_ref <- .joinmefit_association_transform_matrix(
                x,
                term_key,
                term = term,
                x_grid = 0,
                n_draws = length(coeff_draws),
                seed = seed,
                payload = payload
            )
            tf_mat <- tf_mat - matrix(tf_ref[, 1], nrow = nrow(tf_mat), ncol = ncol(tf_mat))
        }

        curve_mat <- if (identical(association_metric, "hazard")) {
            tf_mat * coeff_draws
        } else {
            tf_mat
        }

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

    median_col <- .quantile_name_from_prob(0.5)
    multi_marker <- length(unique(plot_df$marker)) > 1L
    p <- if (multi_marker) {
        ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data[[median_col]], color = .data$marker, fill = .data$marker, group = .data$marker))
    } else {
        ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data[[median_col]]))
    }

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

    p <- if (multi_marker) {
        p + ggplot2::geom_line(linewidth = prediction_style$linewidth %||% 0.8)
    } else {
        p + ggplot2::geom_line(
            color = prediction_style$color %||% "steelblue",
            linewidth = prediction_style$linewidth %||% 0.8
        )
    }

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

    y_lab <- if (identical(association_metric, "hazard")) {
        paste0(term, " contribution")
    } else {
        paste0(term_key, " transform")
    }

    plot_title <- if (identical(association_metric, "hazard")) {
        paste("Association contribution:", term)
    } else {
        paste("Association transform:", term)
    }

    p <- p + ggplot2::labs(
        title = plot_title,
        x = term_key,
        y = y_lab,
        color = if (multi_marker) "marker" else NULL,
        fill = if (multi_marker) "marker" else NULL
    ) + theme_fn()

    y_scale_adjust <- .joinmefit_assoc_transform_limits(plot_df, association_metric)
    if (!is.null(y_scale_adjust)) {
        p <- p + y_scale_adjust
    }
    p
}

.association_plot_markers <- function(x, term_key, marker = NULL) {
    if (term_key %in% c("corr", "vcov")) {
        return("all")
    }
    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    available_markers <- if (!is.null(x$dataLong) && marker_var %in% names(x$dataLong)) {
        unique(as.character(stats::na.omit(x$dataLong[[marker_var]])))
    } else {
        character(0)
    }
    if (is.null(marker)) {
        if (length(available_markers) == 0) "all" else available_markers
    } else {
        marker
    }
}

.joinmefit_transform_specs <- function(x) {
    x$config$transforms_spec %||% x$call$transforms %||% list()
}

.joinmefit_assoc_component_index <- function(term) {
    term <- as.character(term %||% "")[1]
    if (!grepl("\\[\\d+\\]$", term)) {
        return(1L)
    }
    as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", term))
}

.joinmefit_transform_spec_for_term <- function(x, term_key, term = term_key) {
    specs <- .joinmefit_transform_specs(x)
    specs[[term_key]] %||% list(type = "identity")
}

.joinmefit_plot_raw_knot_range <- function(tf_spec) {
    knots <- as.numeric(tf_spec$knots %||% tf_spec$x %||% numeric(0))
    if (!length(knots) || !all(is.finite(knots))) {
        return(NULL)
    }
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

.joinmefit_support_outside_spline_range <- function(support, tf_spec) {
    knots <- as.numeric(tf_spec$knots %||% numeric(0))
    if (length(support) != 2L || length(knots) < 2L) {
        return(FALSE)
    }
    boundary <- c(knots[1], knots[length(knots)])
    if (.transform_uses_expit_input(tf_spec)) {
        support <- .transform_input_for_spec(support, tf_spec)
    }
    any(support < boundary[1] | support > boundary[2])
}

.default_joinmefit_association_grid <- function(x, term_key, marker = NULL, association_range = NULL, n = 200) {
    if (!is.null(association_range)) {
        return(seq(association_range[1], association_range[2], length.out = n))
    }

    if (identical(term_key, "corr")) {
        support <- .joinmefit_association_support_range(x, term_key = term_key, marker = marker)
        if (!is.null(support)) {
            return(seq(support[1], support[2], length.out = n))
        }
        return(seq(-1, 1, length.out = n))
    }
    if (identical(term_key, "vcov")) {
        support <- .joinmefit_association_support_range(x, term_key = term_key, marker = marker)
        if (!is.null(support)) {
            return(seq(support[1], support[2], length.out = n))
        }
    }

    tf_spec <- .joinmefit_transform_specs(x)[[term_key]]
    tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")
    support <- .joinmefit_association_support_range(x, term_key = term_key, marker = marker)
    if (!is.null(support)) {
        if (.is_ispline_transform_type(tf_type) && !is.null(tf_spec$knots)) {
            kr <- .joinmefit_plot_raw_knot_range(tf_spec)
            if (!is.null(kr) && diff(kr) > 0 && .joinmefit_support_outside_spline_range(support, tf_spec)) {
                cli::cli_warn(c(
                    x = "Model-implied support for {.val {term_key}} extends beyond the fitted spline knot range.",
                    i = "Using knot support [{format(signif(kr[1], 4), scientific = FALSE)}, {format(signif(kr[2], 4), scientific = FALSE)}] to avoid unsupported spline extrapolation."
                ))
                return(seq(kr[1], kr[2], length.out = n))
            }
        }
        return(seq(support[1], support[2], length.out = n))
    }

    response_var <- all.vars(x$formulaLong)[1] %||% "y"
    marker_var <- .joinmefit_call_arg_chr(x$call, "marker_var", "marker")
    data_long <- x$dataLong
    if (!is.null(marker) && marker_var %in% names(data_long)) {
        data_long <- data_long[as.character(data_long[[marker_var]]) %in% as.character(marker), , drop = FALSE]
    }
    y_obs <- data_long[[response_var]]
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

    if (term_key %in% c("cv_total", "cv_mean", "cv_marker") && length(y_obs) > 1L) {
        xr <- stats::quantile(y_obs, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
        if (all(is.finite(xr)) && diff(xr) > 0) {
            if (.is_ispline_transform_type(tf_type) && !is.null(tf_spec$knots)) {
                kr <- .joinmefit_plot_raw_knot_range(tf_spec)
                if (!is.null(kr) && diff(kr) > 0 && .joinmefit_support_outside_spline_range(xr, tf_spec)) {
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

    if (!is.null(tf_spec$x)) {
        xr <- .joinmefit_plot_raw_knot_range(list(type = tf_type, x = tf_spec$x))
        if (all(is.finite(xr)) && diff(xr) > 0) {
            return(seq(xr[1], xr[2], length.out = n))
        }
    }
    if (!is.null(tf_spec$knots)) {
        xr <- .joinmefit_plot_raw_knot_range(tf_spec)
        if (all(is.finite(xr)) && diff(xr) > 0) {
            return(seq(xr[1], xr[2], length.out = n))
        }
    }

    if (length(y_obs) > 1L) {
        xr <- stats::quantile(y_obs, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
        if (all(is.finite(xr)) && diff(xr) > 0) return(seq(xr[1], xr[2], length.out = n))
    }

    seq(-1, 1, length.out = n)
}

.joinmefit_assoc_transform_limits <- function(plot_df, association_metric) {
    if (!identical(association_metric, "transform")) {
        return(NULL)
    }

    q_cols <- grep("^q", names(plot_df), value = TRUE)
    vals <- unlist(plot_df[, unique(c("mean", q_cols)), drop = FALSE], use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
        return(NULL)
    }

    if (min(vals) >= -0.02 && max(vals) <= 1.02) {
        return(ggplot2::coord_cartesian(ylim = c(0, 1)))
    }

    NULL
}

.joinmefit_association_marker_weight_draws <- function(x, term_key, marker_level, n_draws, seed = 1, payload = NULL) {
    if (!(term_key %in% c("cv_total", "cs_total", "cv_marker", "cs_marker"))) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    marker_levels <- as.character(x$stan_data$marker_levels %||% unique(stats::na.omit(x$dataLong$marker)) %||% character(0))
    if (length(marker_levels) == 0L) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    marker_idx <- match(as.character(marker_level), marker_levels)
    if (is.na(marker_idx)) {
        return(matrix(1, nrow = n_draws, ncol = 1L))
    }

    weight_draws <- .joinmefit_marker_weight_draws(x, n_draws = n_draws, seed = seed, payload = payload)
    if (is.null(weight_draws) || ncol(weight_draws) < marker_idx) {
        return(matrix(1 / length(marker_levels), nrow = n_draws, ncol = 1L))
    }

    matrix(weight_draws[, marker_idx, drop = TRUE] / length(marker_levels), ncol = 1L)
}

.joinmefit_marker_weight_draws <- function(x, n_draws, seed = 1, payload = NULL) {
    sd <- x$stan_data
    marker_levels <- as.character(sd$marker_levels %||% unique(stats::na.omit(x$dataLong$marker)) %||% character(0))
    n_markers <- length(marker_levels)
    if (n_markers == 0L) {
        return(NULL)
    }

    if (!is.null(payload$marker_weight_draws)) {
        return(as.matrix(payload$marker_weight_draws[, seq_len(min(n_markers, ncol(payload$marker_weight_draws))), drop = FALSE]))
    }

    eff_names <- paste0("marker_weights_eff[", seq_len(n_markers), "]")
    dmat <- tryCatch(
        .get_draws_matrix(x$fit, variables = eff_names, draws = n_draws, seed = seed),
        error = function(e) NULL
    )
    if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
        return(as.matrix(dmat[, eff_names, drop = FALSE]))
    }

    base_names <- paste0("marker_weights[", seq_len(n_markers), "]")
    dmat <- tryCatch(
        .get_draws_matrix(x$fit, variables = base_names, draws = n_draws, seed = seed),
        error = function(e) NULL
    )
    if (!is.null(dmat) && all(base_names %in% colnames(dmat))) {
        return(as.matrix(dmat[, base_names, drop = FALSE]))
    }

    base_weights <- as.numeric(sd$marker_weights %||% rep(1, n_markers))
    if (length(base_weights) < n_markers) {
        base_weights <- c(base_weights, rep(1, n_markers - length(base_weights)))
    }
    matrix(rep(base_weights[seq_len(n_markers)], each = n_draws), nrow = n_draws, byrow = FALSE)
}

.joinmefit_assoc_channel_map <- function(term_key) {
    switch(term_key,
        cv_total = list(
            eff_prefix = "coeff_cv_eff",
            base_coeff = "coeff_cv",
            n_coeff = "n_coeff_cv",
            knots = "knots_cv",
            degree = "spline_degree_cv"
        ),
        cs_total = list(
            eff_prefix = "coeff_cs_eff",
            base_coeff = "coeff_cs",
            n_coeff = "n_coeff_cs",
            knots = "knots_cs",
            degree = "spline_degree_cs"
        ),
        corr = list(
            eff_prefix = "coeff_corr_eff",
            base_coeff = "coeff_corr",
            n_coeff = "n_coeff_corr",
            knots = "knots_corr",
            degree = "spline_degree_corr"
        ),
        vcov = list(
            eff_prefix = "coeff_vcov_eff",
            base_coeff = "coeff_vcov",
            n_coeff = "n_coeff_vcov",
            knots = "knots_vcov",
            degree = "spline_degree_vcov"
        ),
        cv_mean = list(
            eff_prefix = "coeff_cv_mean_eff",
            base_coeff = "coeff_cv_mean",
            n_coeff = "n_coeff_cv_mean",
            knots = "knots_cv_mean",
            degree = "spline_degree_cv_mean"
        ),
        cv_marker = list(
            eff_prefix = "coeff_cv_marker_eff",
            base_coeff = "coeff_cv_marker",
            n_coeff = "n_coeff_cv_marker",
            knots = "knots_cv_marker",
            degree = "spline_degree_cv_marker"
        ),
        cs_mean = list(
            eff_prefix = "coeff_cs_mean_eff",
            base_coeff = "coeff_cs_mean",
            n_coeff = "n_coeff_cs_mean",
            knots = "knots_cs_mean",
            degree = "spline_degree_cs_mean"
        ),
        cs_marker = list(
            eff_prefix = "coeff_cs_marker_eff",
            base_coeff = "coeff_cs_marker",
            n_coeff = "n_coeff_cs_marker",
            knots = "knots_cs_marker",
            degree = "spline_degree_cs_marker"
        ),
        NULL
    )
}

.joinmefit_association_transform_matrix <- function(x, term_key, term = term_key, x_grid, n_draws, seed = 1, payload = NULL) {
    tf_spec <- .joinmefit_transform_spec_for_term(x, term_key, term = term)
    tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")

    if (identical(tf_type, "identity")) {
        return(matrix(rep(as.numeric(x_grid), each = n_draws), nrow = n_draws))
    }

    if (identical(tf_type, "functional")) {
        tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
        vals <- as.numeric(tf_fun(x_grid))
        return(matrix(rep(vals, each = n_draws), nrow = n_draws))
    }

    if (identical(tf_type, "pwlin")) {
        tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
        vals <- as.numeric(tf_fun(x_grid))
        return(matrix(rep(vals, each = n_draws), nrow = n_draws))
    }

    if (.is_ispline_transform_type(tf_type)) {
        return(.joinmefit_ispline_transform_matrix(x, term_key, term = term, x_grid, n_draws = n_draws, seed = seed, payload = payload))
    }

    tf_fun <- .joinme_make_assoc_transform(tf_spec, term_key)
    vals <- as.numeric(tf_fun(x_grid))
    matrix(rep(vals, each = n_draws), nrow = n_draws)
}

.joinmefit_ispline_transform_matrix <- function(x, term_key, term = term_key, x_grid, n_draws, seed = 1, payload = NULL) {
    if (!requireNamespace("splines2", quietly = TRUE)) {
        cli::cli_abort(c(
            x = "Package {.pkg splines2} is required for spline-based association plotting.",
            i = "Install {.pkg splines2} to plot spline transforms."
        ))
    }

    sd <- x$stan_data
    map <- .joinmefit_assoc_channel_map(term_key)
    tf_spec <- .joinmefit_transform_spec_for_term(x, term_key, term = term)

    knots <- as.numeric(sd[[map$knots]] %||% tf_spec$knots %||% tf_spec$x)
    degree <- as.integer(sd[[map$degree]] %||% tf_spec$degree %||% 3L)
    n_coeff <- as.integer(sd[[map$n_coeff]] %||% length(sd[[map$base_coeff]] %||% tf_spec$coeff %||% numeric(0)))
    coeff_draws <- .joinmefit_transform_coeff_draws(
        x = x,
        eff_prefix = map$eff_prefix,
        base_coeff = sd[[map$base_coeff]] %||% tf_spec$coeff,
        n_coeff = n_coeff,
        component_index = if (term_key %in% c("corr", "vcov")) .joinmefit_assoc_component_index(term) else NULL,
        n_draws = n_draws,
        seed = seed,
        payload = payload,
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

.joinmefit_transform_coeff_draws <- function(x, eff_prefix, base_coeff, n_coeff, n_draws, seed = 1, payload = NULL, term_key = NULL, term = term_key, component_index = NULL) {
    n_coeff <- as.integer(n_coeff %||% 0L)
    if (n_coeff < 1L) {
        return(matrix(0, nrow = n_draws, ncol = 0L))
    }

    payload_key <- term %||% term_key
    if (!is.null(payload_key) && !is.null(payload$transform_coeff_draws[[payload_key]])) {
        return(as.matrix(payload$transform_coeff_draws[[payload_key]]))
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

.joinmefit_assoc_coeff_draws <- function(x, term, payload = NULL, seed = 1) {
    term_key <- if (grepl("^corr", term)) {
        "corr"
    } else if (grepl("^vcov", term)) {
        "vcov"
    } else {
        term
    }
    if (!is.null(payload$coeff_draws[[term_key]])) {
        vals <- payload$coeff_draws[[term_key]]
        if (is.matrix(vals)) {
            target_var <- NULL
            if (!is.null(payload$term_map) && nrow(payload$term_map) > 0L) {
                target_rows <- payload$term_map[payload$term_map$term == term, , drop = FALSE]
                if (nrow(target_rows) > 0L) {
                    target_var <- as.character(target_rows$variable[[1L]])
                }
            }
            if (is.null(target_var) && grepl("^(corr|vcov)\\[\\d+\\]$", term)) {
                target_var <- paste0("alpha_", term)
            }
            if (!is.null(target_var) && !is.null(colnames(vals)) && target_var %in% colnames(vals)) {
                return(as.numeric(vals[, target_var, drop = TRUE]))
            }
            return(as.numeric(vals[, 1]))
        }
        return(as.numeric(vals))
    }

    assoc_draws <- extract.JoinMeFit(x, what = "assoc", term = term, keep_chains = FALSE)$draws
    as.numeric(assoc_draws[, 1])
}

.joinmefit_get_association_plot_payload <- function(x, seed = 1) {
    payload <- x$config$association_plot_payload %||% x$cache_get("association_plot_payload")
    if (!is.null(payload)) {
        return(payload)
    }
    if (is.null(x$fit)) {
        return(NULL)
    }

    payload <- tryCatch(
        .build_joinmefit_association_plot_payload(
            fit = x$fit,
            stan_data = x$stan_data,
            config = x$config,
            dataLong = x$dataLong,
            seed = seed
        ),
        error = function(e) NULL
    )
    if (!is.null(payload)) {
        x$config$association_plot_payload <- payload
        x$cache_set("association_plot_payload", payload)
    }
    payload
}

.joinmefit_association_support_range <- function(x, term_key, marker = NULL) {
    payload <- .joinmefit_get_association_plot_payload(x)
    support <- payload$support
    if (is.null(support) || !nrow(support)) {
        return(NULL)
    }

    marker_key <- if (is.null(marker) || (length(marker) == 1L && is.na(marker))) "all" else as.character(marker[[1]])
    rows <- support[support$term == term_key & support$marker %in% c(marker_key, "all"), , drop = FALSE]
    if (!nrow(rows)) {
        return(NULL)
    }
    c(min(rows$lower, na.rm = TRUE), max(rows$upper, na.rm = TRUE))
}

.build_joinmefit_association_plot_payload <- function(fit, stan_data, config, dataLong, seed = 1) {
    assoc_flags <- config$assoc %||% list()
    active_terms <- names(assoc_flags)[vapply(assoc_flags, function(flag) isTRUE(as.logical(flag)), logical(1))]
    if (!length(active_terms)) {
        return(NULL)
    }

    payload <- list(
        coeff_draws = list(),
        marker_weight_draws = NULL,
        transform_coeff_draws = list(),
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
            corr_vars <- grep(paste0("^alpha_", term_key, "\\["), posterior::variables(.get_draws_obj(fit)), value = TRUE)
            if (length(corr_vars)) {
                payload$coeff_draws[[term_key]] <- .get_draws_matrix(fit, variables = corr_vars, seed = seed)[, corr_vars, drop = FALSE]
                payload$term_map <- rbind(
                    payload$term_map,
                    data.frame(
                        term = sub("^alpha_", "", corr_vars),
                        variable = corr_vars,
                        stringsAsFactors = FALSE
                    )
                )
            }
        } else {
            alpha_var <- alpha_map[[term_key]]
            if (!is.null(alpha_var)) {
                payload$coeff_draws[[term_key]] <- .get_draws_matrix(fit, variables = alpha_var, seed = seed)[, 1, drop = FALSE]
                payload$term_map <- rbind(payload$term_map, data.frame(term = term_key, variable = alpha_var, stringsAsFactors = FALSE))
            }
        }

        tf_spec <- payload$transform_specs[[term_key]] %||% list(type = "identity")
        tf_type <- .canonicalise_transform_type(tf_spec$type %||% "identity")
        if (.is_ispline_transform_type(tf_type)) {
            map <- .joinmefit_assoc_channel_map(term_key)
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
                            payload$transform_coeff_draws[[component_terms[[m]]]] <- as.matrix(dmat[, eff_names, drop = FALSE])
                        }
                    }
                } else {
                    eff_names <- paste0(map$eff_prefix, "[", seq_len(n_coeff), "]")
                    dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
                    if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                        payload$transform_coeff_draws[[term_key]] <- as.matrix(dmat[, eff_names, drop = FALSE])
                    }
                }
            }
        }
    }

    if (any(active_terms %in% c("cv_total", "cs_total", "cv_marker", "cs_marker"))) {
        n_markers <- as.integer(stan_data$D %||% length(stan_data$marker_levels %||% numeric(0)))
        if (n_markers > 0L) {
            eff_names <- paste0("marker_weights_eff[", seq_len(n_markers), "]")
            dmat <- tryCatch(.get_draws_matrix(fit, variables = eff_names, seed = seed), error = function(e) NULL)
            if (!is.null(dmat) && all(eff_names %in% colnames(dmat))) {
                payload$marker_weight_draws <- as.matrix(dmat[, eff_names, drop = FALSE])
            }
        }
    }

    payload$support <- .joinmefit_model_implied_support(fit, stan_data, config, dataLong, seed = seed)
    payload
}

.joinmefit_model_implied_support <- function(fit, stan_data, config, dataLong, seed = 1) {
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
    beta_vars <- paste0("beta_used_in_likelihood[", seq_len(P), "]")
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
                    w_vars <- c(w_vars, paste0("w_idscaled[", i, ",", d, ",", q, "]"))
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
        out[[length(out) + 1L]] <- data.frame(term = "corr", marker = "all", lower = -1, upper = 1, source = "theoretical", stringsAsFactors = FALSE)
    }
    if ("vcov" %in% terms && isTRUE((stan_data$Q_idm %||% 0L) > 0L)) {
        vcov_map <- .assoc_cov_feature_map(
            q_idm = as.integer(stan_data$Q_idm),
            diagonal_only = as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L,
            include_diag = TRUE
        )
        alpha_mean <- mean_of_vars(paste0("alpha_L[", seq_len(nrow(vcov_map)), "]"), default = 0)
        lambda_mean <- mean_of_vars(paste0("lambda_L[", seq_len(nrow(vcov_map)), "]"), default = 0)
        tau_l_mean <- mean_of_vars("tau_L", default = 0)[1]
        z_l_mean <- mean_of_vars(paste0("z_L[", seq_len(n_id), "]"), default = 0)
        k_cov <- as.integer(stan_data$K_cov %||% 0L)
        beta_mean <- if (k_cov > 0L && nrow(vcov_map) > 0L) {
            beta_flat <- mean_of_vars(as.vector(outer(seq_len(nrow(vcov_map)), seq_len(k_cov), function(m, k) paste0("beta_L[", m, ",", k, "]"))), default = 0)
            matrix(beta_flat, nrow = nrow(vcov_map), byrow = TRUE)
        } else {
            matrix(0, nrow = nrow(vcov_map), ncol = 0L)
        }
        xcov <- as.matrix(stan_data$Xcov %||% matrix(0, nrow = n_id, ncol = k_cov))
        li_terms <- lapply(seq_len(n_id), function(i) {
            li <- matrix(0, nrow = as.integer(stan_data$Q_idm), ncol = as.integer(stan_data$Q_idm))
            if (nrow(vcov_map) == 0L) {
                return(numeric(0))
            }
            for (m in seq_len(nrow(vcov_map))) {
                lp <- alpha_mean[m] +
                    if (k_cov > 0L) sum(beta_mean[m, ] * xcov[i, ]) else 0 +
                    lambda_mean[m] * tau_l_mean * z_l_mean[min(i, length(z_l_mean))]
                r_ <- vcov_map[m, 1]
                c_ <- vcov_map[m, 2]
                li[r_, c_] <- if (r_ == c_) {
                    if (as.integer(stan_data$corr_diag_link %||% 0L) == 1L) exp(lp) else log1p(exp(lp))
                } else {
                    lp
                }

                    if (identical(association_metric, "hazard") && term_key %in% c("corr", "vcov")) {
                        tf_ref <- .joinmefit_association_transform_matrix(
                            x,
                            term_key,
                            term = term,
                            x_grid = 0,
                            n_draws = length(coeff_draws),
                            seed = seed,
                            payload = payload
                        )
                        tf_mat <- tf_mat - matrix(tf_ref[, 1], nrow = nrow(tf_mat), ncol = ncol(tf_mat))
                    }
            }
            .assoc_vcov_features_from_chol(
                li,
                diagonal_only = as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
            )
        })
        li_vals <- unlist(li_terms, use.names = FALSE)
        li_vals <- li_vals[is.finite(li_vals)]
        if (length(li_vals) > 1L) {
            out[[length(out) + 1L]] <- make_support_rows("vcov", li_vals)
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
        bc <- parse_transform_expr(spec$expr)
        return(function(x) {
            eval_bytecode_vector(
                x = x,
                bytecode = bc$bytecode %||% bc$opcodes,
                const_data = bc$const_data %||% numeric(0)
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
