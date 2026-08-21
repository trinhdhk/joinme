/**
 * @file survival/generated_quantities/dynamic_prediction_declarations.stan
 * @brief Declare conditional survival and cumulative-hazard arrays.
 *
 * @details Both quantities are retained because survival is the exponential of
 * negative cumulative hazard, while diagnostics may require the cumulative
 * hazard on its additive scale.
 */
matrix[n_draws, n_times_surv] surv_prob; // conditional survival probability by draw and future time
matrix[n_draws, n_times_surv] cumhaz_cond; // conditional cumulative hazard by draw and future time

