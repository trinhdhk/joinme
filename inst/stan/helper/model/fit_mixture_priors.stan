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
 * combined class request still has G components,
 * not the Cartesian product of separate G-component mixtures.
 *
 * The component density is placed on the non-centred coordinates themselves.
 * The longitudinal likelihood later sees `L_u * z_u` and `L_v * z_v`, whilst
 * covariance regression sees `alpha_L + beta_L * X + lambda_L .* z_L`.
 * This is the same order used by the R generator. It also explains why the
 * component scales can compensate for the ordinary random-effect scales or
 * covariance loadings. The likelihood may identify their product more readily
 * than the separate factors, producing weak scale identification, a long
 * posterior ridge and low E-BFMI. That geometry is part of the fitted
 * statistical model and is not removed here by hidden anchoring or rescaling.
 * The R-side `diagnosis()` method reports this trade-off explicitly.
 */

if (use_mixture == 1) {
  // STEP 1: assign priors to parameters which describe the common mixture.
  //
  // These statements do not yet evaluate any subject or marker random effect.
  // They describe the baseline class simplex, component centres, component
  // scales, and optional class-membership regression coefficients. The
  // simulator treats these as fixed generative inputs; fitting estimates them.
  // The baseline simplex is common to every selected level. formulaClass
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
  // Every positive component scale receives the original Exponential(1)
  // prior. This statement is intentionally vectorised: it is exactly
  // equivalent to one independent exponential density for every class and
  // packed coordinate, including the normalising constant.
  to_vector(mix_scale) ~ exponential(1);
  mix_class_coefficient_subject ~ normal(0, class_regression_scale);
  mix_class_coefficient_marker ~ normal(0, class_regression_scale);

  // STEP 2: replace the ordinary prior in the subject allocation domain.
  //
  // For a subject i and class g, R generates one label C_i=g, then draws the
  // selected z_u[i] and selected z_L[i] coordinates conditionally independently
  // from component g. Their conditional log densities must therefore be added
  // before marginalising C_i. A second log_sum_exp would incorrectly create
  // independent subject and covariance labels and G x G combined components.
  if (mix_subject == 1 || mix_covariance == 1) {
    for (subject in 1 : n_id) {
      // Start with log Pr(C_i=g | X_i). `latent_progress_probability` uses the
      // same reference-class multinomial-logit calculation as
      // `.mixture_class_probability()` in the R simulator.
      vector[n_classes] component_log_density = log(
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
        // Copy only the requested source coordinates. Their order is the order
        // recorded by R in mix_idx_subject and is paired with the packed
        // segment beginning at mix_start_subject.
        for (coordinate in 1 : mix_dim_subject) {
          selected_subject_effect[coordinate] =
            z_u[subject][mix_idx_subject[coordinate]];
        }
        for (group in 1 : n_classes) {
          // Add log p(z_u[selected] | C_i=group). Multiplication by the
          // established random-effect weight preserves JoiNMe's weighted
          // pseudo-posterior convention; all weights equal one in an ordinary
          // generative analysis such as the worked simulation vignette.
          component_log_density[group] += re_weight_id[subject]
            * latent_progress_component_lpdf(
                selected_subject_effect |
                mix_location[group],
                mix_scale[group],
                mix_start_subject,
                shrinkage
              );
        }
        // `fit_priors.stan` has already added this exact standard-Normal
        // density. Subtract it once so the final target contains the mixture
        // density rather than both the ordinary and mixture priors.
        target += -re_weight_id[subject]
          * std_normal_lpdf(selected_subject_effect);
      }

      if (mix_covariance == 1) {
        vector[mix_dim_covariance] selected_covariance_effect; // Role: selected covariance effect.
        // z_L has one coordinate for every lower-triangular covariance
        // predictor. R constructs the same row-major lower-triangle map and
        // stores its selected source positions in mix_idx_covariance.
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
        // Covariance-regression latent terms also have a standard-Normal
        // ordinary prior. The shrinkage switch changes only their selected
        // replacement component family.
        target += -re_weight_L[subject]
          * std_normal_lpdf(selected_covariance_effect);
      }

      // This is the exact analytic marginalisation
      // log sum_g Pr(C_i=g | X_i) p(z_u,i, z_L,i | C_i=g).
      target += log_sum_exp(component_log_density);
    }
  }

  // STEP 3: replace the ordinary prior in the marker allocation domain.
  //
  // R samples one label per marker and all selected z_v coordinates for that
  // marker use it. This loop is separate from the subject loop because marker
  // and subject allocations are distinct natural units, though they share the
  // same class labels, component rows, and baseline simplex.
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
      );

      vector[mix_dim_marker] selected_marker_effect; // Role: selected marker effect.
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
      target += -re_weight_marker[marker_index]
        * std_normal_lpdf(selected_marker_effect);

      target += log_sum_exp(component_log_density);
    }
  }
}
