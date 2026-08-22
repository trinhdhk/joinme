/**
 * @file include/submodels/latent_class/model/dynamic_prediction.stan
 * @brief Restores the fitted latent-progress distribution during prediction.
 *
 * @details The shared conditioning likelihood supplies the established
 * standard-Normal priors for newly sampled random effects. This module removes
 * those densities on selected coordinates and adds the fitted marginal
 * mixture, preserving one allocation per natural domain.
 */

// -----------------------------------------------------------------------
// Restore the fitted latent-progress prior on dynamically sampled effects.
//
// `partial_draw()` deliberately retains the ordinary standard-Normal
// densities used by established JoiNMe predictions.  The following
// density-ratio adjustment subtracts that density only for selected
// coordinates, then adds the fitted marginal mixture.  This arrangement
// leaves ordinary predictions unchanged and keeps the large threaded
// likelihood helper focused on observation and event calculations.
// -----------------------------------------------------------------------
if (use_dynamic_mixture == 1) {
  if (dynamic_mix_subject == 1 || dynamic_mix_covariance == 1) {
    for (fitted_draw in 1 : n_draws) {
      vector[dynamic_n_classes] component_log_density =
        log(to_vector(
          dynamic_mix_probability_subject[fitted_draw, ]
        ));

      if (dynamic_mix_subject == 1) {
        vector[dynamic_mix_dim_subject] selected_subject_effect; // Role: selected subject effect.
        for (coordinate in 1 : dynamic_mix_dim_subject) {
          selected_subject_effect[coordinate] =
            z_u[fitted_draw][dynamic_mix_idx_subject[coordinate]];
        }
        for (group in 1 : dynamic_n_classes) {
          component_log_density[group] +=
            latent_progress_component_lpdf(
              selected_subject_effect |
              dynamic_mix_location[fitted_draw][group],
              dynamic_mix_scale[fitted_draw][group],
              dynamic_mix_start_subject,
              dynamic_mix_component_family,
              dynamic_mix_component_df
            );
        }
        target += -std_normal_lpdf(selected_subject_effect);
      }

      if (dynamic_mix_covariance == 1) {
        vector[dynamic_mix_dim_covariance] selected_covariance_effect; // Role: selected covariance effect.
        for (coordinate in 1 : dynamic_mix_dim_covariance) {
          selected_covariance_effect[coordinate] =
            z_L[fitted_draw][dynamic_mix_idx_covariance[coordinate]];
        }
        for (group in 1 : dynamic_n_classes) {
          component_log_density[group] +=
            latent_progress_component_lpdf(
              selected_covariance_effect |
              dynamic_mix_location[fitted_draw][group],
              dynamic_mix_scale[fitted_draw][group],
              dynamic_mix_start_covariance,
              dynamic_mix_component_family,
              dynamic_mix_component_df
            );
        }
        target += -std_normal_lpdf(selected_covariance_effect);
      }

      // One log-sum-exp gives the new subject one allocation shared by
      // every requested subject-indexed block.
      target += log_sum_exp(component_log_density);
    }
  }

  if (dynamic_mix_marker == 1) {
    for (fitted_draw in 1 : n_draws) {
      for (marker_index in 1 : n_marker_types) {
        vector[dynamic_n_classes] component_log_density =
          log(to_vector(
            dynamic_mix_probability_marker[fitted_draw][marker_index, ]
          ));
        vector[dynamic_mix_dim_marker] selected_marker_effect; // Role: selected marker effect.

        for (coordinate in 1 : dynamic_mix_dim_marker) {
          selected_marker_effect[coordinate] =
            z_v[fitted_draw, marker_index][
              dynamic_mix_idx_marker[coordinate]
            ];
        }
        for (group in 1 : dynamic_n_classes) {
          component_log_density[group] +=
            latent_progress_component_lpdf(
              selected_marker_effect |
              dynamic_mix_location[fitted_draw][group],
              dynamic_mix_scale[fitted_draw][group],
              dynamic_mix_start_marker,
              dynamic_mix_component_family,
              dynamic_mix_component_df
            );
        }

        target += -std_normal_lpdf(selected_marker_effect);
        // This log-sum-exp gives the marker one allocation.
        target += log_sum_exp(component_log_density);
      }
    }
  }
}
