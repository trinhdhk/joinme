/**
 * @file marker_weight/transformed_parameters/fit.stan
 * @brief Construct effective marker weights.
 *
 * @details
 * Each effective weight is the supplied offset plus its set-level fitted
 * common location and a marker-specific centred unit-scale departure. No
 * additional departure scale is estimated because it would be confounded
 * with the association slope that multiplies the weighted marker aggregate.
 * All Student-t weight sets additionally share one fitted degrees-of-freedom
 * value above two when the R family declaration requests moving df. Shared
 * terms point to the same weight set. Fixed weights have no
 * fitted locations, departures, or degrees of freedom.
 */
  /* -------------------- marker weights (signed and partially pooled) */
  // Goal:
  // - Allow positive and negative marker contributions.
  // - Support either one shared marker-weight structure across all weighted
  //   marker-based association terms or one structure per active term.
  // - Estimate one common location per weight set. This location is informed
  //   jointly by every marker in its set and therefore plays the same role as
  //   a fixed intercept in a hierarchical regression.
  // - Keep marker-specific departures directly interpretable by adding signed
  //   standardised deviations to the supplied offset and fitted common location.
  vector[n_marker_weight_sets] marker_weight_mean = rep_vector(0.0, n_marker_weight_sets); // fitted common location of each set, zero for fixed or inactive weights
  vector[prior_marker_weight_family == 1 && estimate_marker_weight_df == 1 && n_marker_weight_means > 0 ? 1 : 0] marker_weight_df = marker_weight_df_excess + rep_vector(2.0, prior_marker_weight_family == 1 && estimate_marker_weight_df == 1 && n_marker_weight_means > 0 ? 1 : 0); // single learned Student-t degrees of freedom shared by all active weight sets; absent when df is fixed or another family is used
  matrix[n_marker_weight_sets, D] z_marker_weight_sets = rep_matrix(0.0, n_marker_weight_sets, D); // realised marker-specific departures on the fixed unit-scale family coordinate
  vector[D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets] marker_weight_prior_effect = joinme_prior_transform(
    z_marker_weights,
    prior_marker_weight_family,
    rep_vector(0, D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets),
    rep_vector(1, D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets),
    horseshoe_local_marker_weight,
    horseshoe_global_marker_weight,
    horseshoe_slab_marker_weight,
    prior_marker_weight_slab_scale
  ); // marker-specific departures after the selected centred unit-scale family transformation; these enter effective weights directly
  vector[D] marker_weights_eff_cv_total = marker_weights_cv_total; // Role: marker weights eff current value total.
  vector[D] marker_weights_eff_cs_total = marker_weights_cs_total; // Role: marker weights eff current slope total.
  vector[D] marker_weights_eff_cv_marker = marker_weights_cv_marker; // Role: marker weights eff current value marker.
  vector[D] marker_weights_eff_cs_marker = marker_weights_cs_marker; // Role: marker weights eff current slope marker.
  {
    // First form common locations and marker departures before combining them
    // with the supplied offsets into effective association weights.
    if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1) {
      for (s in 1:n_marker_weight_sets) {
        int start_pos = (s - 1) * D + 1; // Role: starting pos.
        int end_pos = s * D; // Role: ending pos.
        int alpha_pos = 6 + M_corr + M_vcov + s; // packed association-block position of this set's independently prior-governed common location
        marker_weight_mean[s] = alpha_prior_effect[alpha_pos];
        z_marker_weight_sets[s] = to_row_vector(marker_weight_prior_effect[start_pos:end_pos]); // direct unit-scale departure z_sd used by every effective weight in set s; alpha supplies association magnitude
      }
    }

    if (assoc_cv_total == 1 && marker_weight_set_cv_total > 0) {
      marker_weights_eff_cv_total = marker_weights_cv_total + marker_weight_mean[marker_weight_set_cv_total] + to_vector(z_marker_weight_sets[marker_weight_set_cv_total]);
    }
    if (assoc_cs_total == 1 && marker_weight_set_cs_total > 0) {
      marker_weights_eff_cs_total = marker_weights_cs_total + marker_weight_mean[marker_weight_set_cs_total] + to_vector(z_marker_weight_sets[marker_weight_set_cs_total]);
    }
    if (assoc_cv_marker == 1 && marker_weight_set_cv_marker > 0) {
      marker_weights_eff_cv_marker = marker_weights_cv_marker + marker_weight_mean[marker_weight_set_cv_marker] + to_vector(z_marker_weight_sets[marker_weight_set_cv_marker]);
    }
    if (assoc_cs_marker == 1 && marker_weight_set_cs_marker > 0) {
      marker_weights_eff_cs_marker = marker_weights_cs_marker + marker_weight_mean[marker_weight_set_cs_marker] + to_vector(z_marker_weight_sets[marker_weight_set_cs_marker]);
    }
  }
