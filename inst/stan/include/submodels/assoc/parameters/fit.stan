/**
 * @file assoc/parameters/fit.stan
 * @brief Parameters for longitudinal--event association coefficients.
 *
 * @details
 * The six scalar channels and the correlation and covariance channels are sampled on standardised scales.  The full-model coefficient programme transforms them under the association slope prior.
 */
  /* association coefficients (non-centred for totals/mean/marker) */
  real z_alpha_cv_total;               // latent for total CV association
  real z_alpha_cs_total;               // latent for total CS association
  real z_alpha_cv_mean;               // latent for mean CV association
  real z_alpha_cs_mean;               // latent for mean CS association
  real z_alpha_cv_marker;              // latent for marker CV association
  real z_alpha_cs_marker;              // latent for marker CS association
  vector[M_corr] z_alpha_corr;        // standardized corr association coefficient latents
  vector[M_vcov] z_alpha_vcov;        // standardized vcov association coefficient latents




