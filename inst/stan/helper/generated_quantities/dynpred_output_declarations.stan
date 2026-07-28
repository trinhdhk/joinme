/**
 * @file dynpred_output_declarations.stan
 * @brief Common dynamic-prediction output declarations.
 *
 * @details These arrays are produced for both ordinary and mixture models.
 * Calculations are separated so mixture-only declarations can be placed before
 * any generated-quantity statement.
 */

/* Outputs per draw and per row/time */
matrix[n_draws, n_obs_long] y_fit_linpred; // Role: outcome fit linpred.
matrix[n_draws, n_obs_long] y_fit_epred; // Role: outcome fit epred.
matrix[n_draws, n_obs_pred] y_pred_linpred; // Role: outcome predictions linpred.
matrix[n_draws, n_obs_pred] y_pred_epred; // Role: outcome predictions epred.
matrix[n_draws, n_obs_pred] y_pred; // Role: outcome predictions.
matrix[n_draws, n_times_surv] surv_prob; // Role: survival prob.
matrix[n_draws, n_times_surv] cumhaz_cond; // Role: cumhaz conditioning.
