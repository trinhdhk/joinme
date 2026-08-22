/**
 * @file include/etc/model/regression_priors.stan
 * @brief Assign the common coefficient priors shared across submodels.
 *
 * @details
 * Each packed raw coefficient receives its declared standard density.  Regularised-horseshoe local scales are coefficient-specific, whereas global and slab scales are shared only within one scientific component-role group.
 */
  // -------------------- Priors
  /**
   * Ordinary regression coefficients.
   *
   * Every raw coefficient receives the standard density implied by its own
   * scientific component and intercept/slope role.  The transformed-parameter
   * block supplies its stated location and scale.  A horseshoe adds one local
   * scale per coefficient and one global/slab pair per component-role group,
   * so regularisation is shared only where the prior declaration says it should be.
   */
  target += joinme_mixed_standard_prior_lpdf(
    regression_prior_raw |
    prior_regression_family,
    prior_regression_df
  );
  for (local_scale in 1 : n_regression_horseshoe_local)
    horseshoe_local_regression[local_scale] ~ student_t(
      prior_regression_horseshoe_local_df[local_scale], 0, 1
    );
  for (prior_group in 1 : n_regression_horseshoe_group) {
    horseshoe_global_regression[prior_group] ~ student_t(
      prior_regression_horseshoe_global_df[prior_group],
      0,
      prior_regression_horseshoe_global_scale[prior_group]
    );
    horseshoe_slab_regression[prior_group] ~ inv_gamma(
      0.5 * prior_regression_horseshoe_slab_df[prior_group],
      0.5 * prior_regression_horseshoe_slab_df[prior_group]
    );
  }
