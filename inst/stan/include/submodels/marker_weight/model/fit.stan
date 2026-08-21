/**
 * @file marker_weight/model/fit.stan
 * @brief Assign the marker-weight departure prior.
 *
 * @details
 * The raw marker coordinates receive the selected centred unit-scale family
 * and enter the effective weights directly. There is deliberately no fitted
 * marker-weight scale: the association slope already multiplies the weighted
 * marker aggregate, so a second multiplier would not be separately identified.
 * Student-t departures either use an explicitly declared fixed degrees of
 * freedom or learn one value per set under a shifted-Gamma prior. When the regularised
 * horseshoe is selected, its local, global and finite-slab hierarchy is
 * assigned here. Common weight locations are governed by the marker-weight
 * intercept role in the common coefficient programme.
 */
  /**
   * Independent standardised prior family for marker-weight departures.
   * Raw quantities are zero-centred and unit-scale; their optional horseshoe
   * transformation is applied in transformed parameters. Their scientific
   * location is then supplied by the declared offset and fitted common mean.
   * This keeps the departure family independent of the class
   * component family, which is declared only by the mixture entry point.
   */
  // Student-t tails use either the explicitly supplied fixed degrees
  // of freedom or one learned value per active set. In the latter case the
  // positive excess nu_s - 2 follows Gamma(shape = 2, rate = 0.1). Adding two
  // gives every conditional Student-t law finite variance while allowing the
  // data from all markers in the set to inform its common tail thickness.
  if (prior_marker_weight_family == 1 && n_marker_weight_means > 0) {
    if (estimate_marker_weight_df == 1) {
      marker_weight_df_excess ~ gamma(2, 0.1);
      for (s in 1:n_marker_weight_means) {
        int start_pos = (s - 1) * D + 1; // first marker departure governed by the learned set-specific degrees of freedom
        int end_pos = s * D; // final marker departure governed by the learned set-specific degrees of freedom
        z_marker_weights[start_pos:end_pos] ~ student_t(marker_weight_df[s], 0, 1);
      }
    } else {
      z_marker_weights ~ student_t(prior_marker_weight_df, 0, 1); // all sets use the degrees of freedom explicitly declared by the analyst
    }
  } else {
    target += joinme_standard_prior_lpdf(
      z_marker_weights | prior_marker_weight_family, prior_marker_weight_df
    );
  }

  // The horseshoe family retains its existing local, global and slab
  // hierarchy and does not instantiate a Student-t degrees-of-freedom model.
  if (prior_marker_weight_family == 4) {
    horseshoe_local_marker_weight ~ student_t(prior_marker_weight_df, 0, 1);
    horseshoe_global_marker_weight ~ student_t(prior_marker_weight_global_df, 0, prior_marker_weight_global_scale);
    horseshoe_slab_marker_weight ~ inv_gamma(0.5 * prior_marker_weight_slab_df, 0.5 * prior_marker_weight_slab_df);
  }
