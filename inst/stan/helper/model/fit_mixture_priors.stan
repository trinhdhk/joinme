/**
 * @file fit_mixture_priors.stan
 * @brief Replace selected standardised random-effect priors by a finite mixture.
 *
 * @details
 * `fit_priors.stan` first contributes the established JoiNMe priors.  This
 * block then subtracts the density of each selected coordinate and adds the
 * requested finite-mixture density.  The resulting target is exactly the
 * mixture model; it is not an additional penalty layered on top of the
 * ordinary prior.
 *
 * Subject and covariance-regression coordinates share one class allocation
 * per subject. Marker effects have one allocation per marker. Consequently a
 * combined clustering request still has G components,
 * not the Cartesian product of separate G-component mixtures.
 */

if (use_mixture == 1) {
  // The baseline simplex is common to every selected level. formulaCluster
  // covariates modify it within the relevant allocation domain.
  if (mix_ordering == 2) {
    // The additive-log-ratio Jacobian makes the Dirichlet statement a density
    // on the ordered simplex rather than on its logit coordinates.
    target += dirichlet_lpdf(mix_probability | mix_probability_prior)
      + sum(log(mix_probability));
  } else {
    mix_probability_unordered[1] ~ dirichlet(mix_probability_prior);
  }
  if (mix_ordered_location_coordinate > 0) {
    mix_location_ordered[1] ~ std_normal();
  }
  to_vector(mix_location_unordered) ~ std_normal();
  to_vector(mix_scale) ~ exponential(1);
  mix_class_coefficient_subject ~ normal(0, class_regression_scale);
  mix_class_coefficient_marker ~ normal(0, class_regression_scale);

  // Subject-indexed domain: one log-sum-exp combines all requested
  // subject and covariance-regression coordinates for a given subject.
  if (mix_subject == 1 || mix_covariance == 1) {
    for (subject in 1 : n_id) {
      vector[n_clusters] component_log_density = log(
        latent_progress_probability(
          mix_probability,
          X_class_subject[subject],
          mix_class_coefficient_subject,
          class_term_start_subject,
          class_term_count_subject
        )
      );

      if (mix_subject == 1) {
        vector[mix_dim_subject] selected_subject_effect; // Role: selected subject effect.
        for (coordinate in 1 : mix_dim_subject) {
          selected_subject_effect[coordinate] =
            z_u[subject][mix_idx_subject[coordinate]];
        }
        for (group in 1 : n_clusters) {
          component_log_density[group] += re_weight_id[subject]
            * latent_progress_component_lpdf(
                selected_subject_effect |
                mix_location[group],
                mix_scale[group],
                mix_start_subject,
                shrinkage
              );
        }
        // The ordinary subject latent prior is always standard Normal.
        target += -re_weight_id[subject]
          * std_normal_lpdf(selected_subject_effect);
      }

      if (mix_covariance == 1) {
        vector[mix_dim_covariance] selected_covariance_effect; // Role: selected covariance effect.
        for (coordinate in 1 : mix_dim_covariance) {
          selected_covariance_effect[coordinate] =
            z_L[subject][mix_idx_covariance[coordinate]];
        }
        for (group in 1 : n_clusters) {
          component_log_density[group] += re_weight_L[subject]
            * latent_progress_component_lpdf(
                selected_covariance_effect |
                mix_location[group],
                mix_scale[group],
                mix_start_covariance,
                shrinkage
              );
        }
        // Covariance-regression latent terms also have a standard Normal
        // ordinary prior, independently of the marker-weight shrinkage choice.
        target += -re_weight_L[subject]
          * std_normal_lpdf(selected_covariance_effect);
      }

      target += log_sum_exp(component_log_density);
    }
  }

  // Marker-indexed domain: selected marker effects share one marker class.
  if (mix_marker == 1) {
    for (marker_index in 1 : D) {
      vector[n_clusters] component_log_density = log(
        latent_progress_probability(
          mix_probability,
          X_class_marker[marker_index],
          mix_class_coefficient_marker,
          class_term_start_marker,
          class_term_count_marker
        )
      );

      vector[mix_dim_marker] selected_marker_effect; // Role: selected marker effect.
      for (coordinate in 1 : mix_dim_marker) {
        selected_marker_effect[coordinate] =
          z_v[marker_index][mix_idx_marker[coordinate]];
      }
      for (group in 1 : n_clusters) {
        component_log_density[group] += re_weight_marker[marker_index]
          * latent_progress_component_lpdf(
              selected_marker_effect |
              mix_location[group],
              mix_scale[group],
              mix_start_marker,
              shrinkage
            );
      }
      target += -re_weight_marker[marker_index]
        * std_normal_lpdf(selected_marker_effect);

      target += log_sum_exp(component_log_density);
    }
  }
}
