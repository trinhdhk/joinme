/**
 * @file include/etc/data/fitted_random_effect_draw_data.stan
 * @brief Supply realised random effects for prediction of fitted subjects.
 *
 * @details
 * These arrays have the same retained-posterior-draw axis as the population
 * parameters.  When `reuse_fitted_re` is one, prediction evaluates the fitted
 * subject's posterior effects directly and does not estimate another set from
 * the supplied longitudinal history.  Neutral zero and identity arrays are
 * supplied when the ordinary new-subject calculation is requested, which
 * keeps one explicit data contract for both prediction routes.
 */
int<lower=0, upper=1> reuse_fitted_re; // one selects posterior random effects belonging to a subject seen during fitting
array[n_draws] vector[n_random_id] fitted_u_id; // realised subject-level effects paired with every retained fit draw
array[n_draws, n_marker_types] vector[n_random_marker] fitted_v_marker; // realised marker-level effects, including the fitted marker prior transformation
array[n_draws, n_marker_types] vector[n_random_marker_id] fitted_z_w; // marker-by-subject latent effects after any marker cross-correlation
array[n_draws] matrix[n_random_marker_id, n_random_marker_id] fitted_L_i; // fitted subject-specific covariance factor before internal-time row scaling
