/**
 * @file include/submodels/latent_class/generated_quantities/dynamic_prediction_calculations.stan
 * @brief Posterior allocation calculation after dynamic conditioning.
 *
 * @details For every retained parent draw, the code combines the fitted class
 * probability with the selected dynamically sampled coordinates, then
 * normalises the component log densities with `softmax`.
 */

for (fitted_draw in 1 : n_draws) {
  posterior_class_probability_new_marker[fitted_draw] =
    rep_matrix(
      1.0 / dynamic_n_classes,
      n_marker_types,
      dynamic_n_classes
    );
}

if (
  use_dynamic_mixture == 1 &&
  (dynamic_mix_subject == 1 || dynamic_mix_covariance == 1)
) {
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
    }

    posterior_class_probability_new_subject[fitted_draw] =
      to_row_vector(softmax(component_log_density));
  }
}

if (
  use_dynamic_mixture == 1 &&
  dynamic_mix_marker == 1
) {
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

      posterior_class_probability_new_marker[fitted_draw, marker_index] =
        to_row_vector(softmax(component_log_density));
    }
  }
}
