/**
 * @file fit_mixture_outputs.stan
 * @brief Posterior allocations and class-centre random effects.
 *
 * @details
 * These quantities do not sample a discrete class.  They report the exact
 * conditional class probabilities implied by each posterior draw.  The
 * deterministic maximum-probability label is supplied only as a convenient
 * summary; probability vectors should be preferred whenever classification
 * uncertainty matters.
 *
 * The calculation deliberately repeats the class-specific part of the model
 * density. It starts from the membership probability implied by the baseline
 * simplex and `formulaClass`, then adds the conditional density of every
 * selected latent coordinate. The ordinary standard-Normal density which is
 * subtracted in the model block is common to all classes. It therefore
 * cancels from Bayes' rule and must not be subtracted here.
 *
 * Subject and covariance-regression contributions enter the same accumulator
 * before normalisation because the R generator assigns them one shared
 * subject class. Marker effects use their separate marker-indexed allocation.
 * Thus these reported probabilities correspond exactly to the generative
 * labels recorded by `simulate_joinme_mix()`.
 */

matrix[n_id, n_classes] posterior_class_probability_subject =
  rep_matrix(1.0 / n_classes, n_id, n_classes); // posterior Pr(C_i = g | all fitted quantities) for every subject and class
matrix[D, n_classes] posterior_class_probability_marker =
  rep_matrix(1.0 / n_classes, D, n_classes); // posterior Pr(C_d = g | all fitted quantities) for every marker and class
array[n_id] int<lower=1, upper=n_classes> posterior_class_subject =
  rep_array(1, n_id); // maximum-probability subject class used only as a concise label
array[D] int<lower=1, upper=n_classes> posterior_class_marker =
  rep_array(1, D); // maximum-probability marker class used only as a concise label

// Class locations are also expanded back to the fitted block dimensions.  A
// zero denotes a coordinate which retained its ordinary centred prior.
matrix[n_classes, R_id] class_mean_subject = rep_matrix(0.0, n_classes, R_id); // Role: class mean subject.
matrix[n_classes, R_mk] class_mean_marker = rep_matrix(0.0, n_classes, R_mk); // Role: class mean marker.
matrix[n_classes, M_cov] class_mean_covariance = rep_matrix(0.0, n_classes, M_cov); // Role: class mean covariance regression.

if (use_mixture == 1) {
  for (group in 1 : n_classes) {
    if (mix_subject == 1) {
      for (coordinate in 1 : mix_dim_subject) {
        class_mean_subject[group, mix_idx_subject[coordinate]] =
          mix_location[
            group,
            mix_start_subject + coordinate - 1
          ];
      }
    }
    if (mix_marker == 1) {
      for (coordinate in 1 : mix_dim_marker) {
        class_mean_marker[group, mix_idx_marker[coordinate]] =
          mix_location[
            group,
            mix_start_marker + coordinate - 1
          ];
      }
    }
    if (mix_covariance == 1) {
      for (coordinate in 1 : mix_dim_covariance) {
        class_mean_covariance[group, mix_idx_covariance[coordinate]] =
          mix_location[
            group,
            mix_start_covariance + coordinate - 1
          ];
      }
    }
  }

  if (mix_subject == 1 || mix_covariance == 1) {
    for (subject in 1 : n_id) {
      vector[n_classes] component_log_density = log(
        latent_progress_probability(
          mix_probability,
          X_class_subject[subject],
          mix_class_coefficient_subject,
          class_term_start_subject,
          class_term_count_subject
        )
      ); // log membership probability plus the selected latent densities, one entry per candidate subject class

      if (mix_subject == 1) {
        vector[mix_dim_subject] selected_subject_effect; // selected z_u coordinates in the same source order used by simulation and fitting
        for (coordinate in 1 : mix_dim_subject) {
          selected_subject_effect[coordinate] =
            z_u[subject][mix_idx_subject[coordinate]];
        }
        for (group in 1 : n_classes) {
          component_log_density[group] += re_weight_id[subject]
            * latent_progress_component_lpdf(
                selected_subject_effect |
                mix_location[group],
                mix_scale[group],
                mix_start_subject,
                shrinkage
              );
        }
      }

      if (mix_covariance == 1) {
        vector[mix_dim_covariance] selected_covariance_effect; // selected z_L coordinates sharing this subject's candidate class
        for (coordinate in 1 : mix_dim_covariance) {
          selected_covariance_effect[coordinate] =
            z_L[subject][mix_idx_covariance[coordinate]];
        }
        for (group in 1 : n_classes) {
          component_log_density[group] += re_weight_L[subject]
            * latent_progress_component_lpdf(
                selected_covariance_effect |
                mix_location[group],
                mix_scale[group],
                mix_start_covariance,
                shrinkage
              );
        }
      }

      posterior_class_probability_subject[subject] =
        to_row_vector(softmax(component_log_density));
      {
        int most_probable_group = 1; // running index of the subject class with greatest conditional log probability
        for (group in 2 : n_classes) {
          if (
            component_log_density[group] >
              component_log_density[most_probable_group]
          ) {
            most_probable_group = group;
          }
        }
        posterior_class_subject[subject] = most_probable_group;
      }
    }
  }

  if (mix_marker == 1) {
    for (marker_index in 1 : D) {
      vector[n_classes] component_log_density = log(
        latent_progress_probability(
          mix_probability,
          X_class_marker[marker_index],
          mix_class_coefficient_marker,
          class_term_start_marker,
          class_term_count_marker
        )
      ); // log membership probability plus the selected latent densities, one entry per candidate marker class

      vector[mix_dim_marker] selected_marker_effect; // selected z_v coordinates in the same source order used by the R generator
      for (coordinate in 1 : mix_dim_marker) {
        selected_marker_effect[coordinate] =
          z_v[marker_index][mix_idx_marker[coordinate]];
      }
      for (group in 1 : n_classes) {
        component_log_density[group] += re_weight_marker[marker_index]
          * latent_progress_component_lpdf(
              selected_marker_effect |
              mix_location[group],
              mix_scale[group],
              mix_start_marker,
              shrinkage
            );
      }

      posterior_class_probability_marker[marker_index] =
        to_row_vector(softmax(component_log_density));
      {
        int most_probable_group = 1; // running index of the marker class with greatest conditional log probability
        for (group in 2 : n_classes) {
          if (
            component_log_density[group] >
              component_log_density[most_probable_group]
          ) {
            most_probable_group = group;
          }
        }
        posterior_class_marker[marker_index] = most_probable_group;
      }
    }
  }
}
