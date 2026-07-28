/**
 * @file dynpred_mixture_output_declarations.stan
 * @brief Conditional latent-class probability outputs for new units.
 *
 * @details Subject-domain and marker-domain probabilities are retained
 * separately because their allocation units differ, even though they share the
 * same fitted class labels and component probabilities.
 */

/**
 * Conditional class probabilities for the dynamically sampled effects.
 *
 * Rows correspond to retained parent-model draws.  These probabilities are
 * conditional on the newly observed longitudinal history because that history
 * informs `z_u`, `z_L`, and (under the established dynamic programme) `z_v`.
 * Inactive allocation domains retain the neutral probability 1/G.
 */
matrix[n_draws, dynamic_n_clusters]
  posterior_class_probability_new_subject =
    rep_matrix(
      1.0 / dynamic_n_clusters,
      n_draws,
      dynamic_n_clusters
    );
array[n_draws] matrix[n_marker_types, dynamic_n_clusters]
  posterior_class_probability_new_marker;
