/**
 * @file marker_weight/data/dynamic_prediction.stan
 * @brief Declare base and posterior marker weights for dynamic prediction.
 *
 * @details
 * Offset values provide a stable reference, while one matrix per weighted
 * channel supplies the effective marker weights retained from every posterior
 * draw. Those retained values already contain the fitted common location and
 * direct unit-scale departure. Dynamic prediction therefore uses exactly the
 * effective weight retained during fitting, with no further transformation.
 */
  vector[n_marker_types] marker_weights_cv_total;            // reference weights for cv_total
  vector[n_marker_types] marker_weights_cs_total;            // reference weights for cs_total
  vector[n_marker_types] marker_weights_cv_marker;           // reference weights for cv_marker
  vector[n_marker_types] marker_weights_cs_marker;           // reference weights for cs_marker
  matrix[n_draws, n_marker_types] marker_weights_draws_cv_total; // per-draw weights for cv_total
  matrix[n_draws, n_marker_types] marker_weights_draws_cs_total; // per-draw weights for cs_total
  matrix[n_draws, n_marker_types] marker_weights_draws_cv_marker; // per-draw weights for cv_marker
  matrix[n_draws, n_marker_types] marker_weights_draws_cs_marker; // per-draw weights for cs_marker
