/**
 * @file marker_weight/parameters/fit.stan
 * @brief Parameters for common marker-weight locations and marker departures.
 *
 * @details
 * One raw common location is fitted per active weight set. Marker-specific
 * coordinates have a family-specific centred unit-scale law and enter the
 * effective weights without another fitted multiplier. This restriction
 * prevents a marker-weight scale from competing with the association slope
 * for the same multiplicative information. When the family
 * name requests moving Student-t degrees of freedom, each active set learns
 * its own positive excess above two. Optional regularised-horseshoe auxiliaries are confined to
 * this submodel.
 */
  /* common marker-weight locations and marker-specific departures (global across subjects) */
  vector[n_marker_weight_means] z_marker_weight_mean; // standardised common locations, one per estimated marker-weight set
  vector[D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets] z_marker_weights; // latent signed coordinates governed directly by the centred unit-scale marker-weight family
  vector<lower=0>[prior_marker_weight_family == 1 && estimate_marker_weight_df == 1 && n_marker_weight_means > 0 ? n_marker_weight_means : 0] marker_weight_df_excess; // learned nu minus two for each active moving-df Student-t weight set; positivity guarantees finite conditional variance
  vector<lower=0>[prior_marker_weight_family == 4 ? D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets : 0] horseshoe_local_marker_weight; // local scales for marker weights
  vector<lower=0>[prior_marker_weight_family == 4 ? 1 : 0] horseshoe_global_marker_weight; // global scale for marker weights
  vector<lower=0>[prior_marker_weight_family == 4 ? 1 : 0] horseshoe_slab_marker_weight; // slab auxiliary for marker weights
