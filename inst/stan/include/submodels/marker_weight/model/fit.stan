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
 * freedom or learn one value shared by all sets under a shifted-Gamma prior. When the regularised
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
  // of freedom or one learned value shared by every active set. In the latter
  // case the positive excess nu - 2 follows Gamma(shape = 2, rate = 0.1). Adding two
  // gives every conditional Student-t law finite variance while allowing the
  // departures from every marker and every weight set to inform one common
  // tail thickness.
  if (prior_marker_weight_family == 1 && n_marker_weight_means > 0) {
    if (estimate_marker_weight_df == 1) {
      marker_weight_df_excess ~ gamma(2, 0.1);
      z_marker_weights ~ student_t(marker_weight_df[1], 0, 1); // every marker departure in every active set shares the single learned tail parameter
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
