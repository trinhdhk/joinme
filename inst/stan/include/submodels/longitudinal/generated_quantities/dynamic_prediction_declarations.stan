/**
 * @file longitudinal/generated_quantities/dynamic_prediction_declarations.stan
 * @brief Declare fitted and future longitudinal prediction arrays.
 *
 * @details Linear predictors, expected responses and posterior predictive
 * responses are kept distinct so plotting and validation can select the
 * appropriate longitudinal estimand.
 */
matrix[n_draws, n_obs_long] y_fit_linpred; // longitudinal linear predictors at observed rows
matrix[n_draws, n_obs_long] y_fit_epred; // expected observed longitudinal responses
matrix[n_draws, n_obs_pred] y_pred_linpred; // longitudinal linear predictors at future rows
matrix[n_draws, n_obs_pred] y_pred_epred; // expected future longitudinal responses
matrix[n_draws, n_obs_pred] y_pred; // posterior predictive future longitudinal responses

