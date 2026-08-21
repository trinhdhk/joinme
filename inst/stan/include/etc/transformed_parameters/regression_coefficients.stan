/**
 * @file include/etc/transformed_parameters/regression_coefficients.stan
 * @brief Pack and transform coefficients shared across scientific submodels.
 *
 * @details
 * This first transformed-parameter fragment reconstructs the common coefficient vector from raw parameters, applies coefficient-specific prior transformations once, and restores longitudinal, survival, association, functional and covariance-regression names.  It must precede every component-specific transformed-parameter fragment.
 */
  /**
  * @brief Time-scaled coefficient vectors and derived random effects for both longitudinal and survival models.
   *
   * This block constructs:
   * 1. Scaled versions of fixed and random effect SDs where time transformations are applied
   * 2. Cholesky factors for random effect covariance structures
   * 3. Marker-average effects used in survival model
  * 4. ID-specific variance structure for covariance regression
  */

  /**
   * STEP 1: collect every ordinary regression coefficient into the shared
   * prior layout prepared by R.
   *
   * The fitted model retains familiar scientific names such as `beta`,
   * `gamma_w` and `beta_sigma`.  Their sampled counterparts are standardised
   * raw parameters.  Packing the raw parameters once permits intercept and
   * slope priors to differ within a component, and permits different
   * components to use different distribution families, without duplicating
   * the same transformation in many places.
   */
  vector[n_regression_prior] regression_prior_raw; // all standardised regression parameters in the R-defined common order
  vector[n_regression_prior] regression_prior_effect; // corresponding coefficients on their scientific model scales
  vector[P] beta; // longitudinal population coefficients on the original time scale
  array[K_event] vector[p_w] gamma_w; // event-covariate slopes on the log-hazard scale
  vector[P_sigma] beta_sigma; // sigma-regression coefficients on its linear-predictor scale
  vector[P_nu] beta_nu; // degrees-of-freedom-regression coefficients on its linear-predictor scale
  vector[P_phi] beta_phi; // dispersion-regression coefficients on its linear-predictor scale
  vector[P_alpha] beta_alpha; // skewness-regression coefficients on its linear-predictor scale
  vector[P_kappa] beta_kappa; // beta-sample-size-regression coefficients on its linear-predictor scale
  vector[P_tau] beta_tau; // quantile-regression coefficients on its linear-predictor scale
  vector[n_alpha_prior] alpha_prior_effect; // association slopes followed by fitted common marker-weight locations

  if (P > 0)
    regression_prior_raw[prior_start_beta : prior_start_beta + P - 1] = z_beta;

  // Association channels have a fixed documented order.  Common marker-weight
  // locations follow them, so these locations can use the marker-weight
  // mean prior while association coefficients use the association slope
  // prior even though both are consumed together downstream.
  regression_prior_raw[prior_start_alpha] = z_alpha_cv_total;
  regression_prior_raw[prior_start_alpha + 1] = z_alpha_cs_total;
  regression_prior_raw[prior_start_alpha + 2] = z_alpha_cv_mean;
  regression_prior_raw[prior_start_alpha + 3] = z_alpha_cs_mean;
  regression_prior_raw[prior_start_alpha + 4] = z_alpha_cv_marker;
  regression_prior_raw[prior_start_alpha + 5] = z_alpha_cs_marker;
  if (M_corr > 0)
    regression_prior_raw[prior_start_alpha + 6 : prior_start_alpha + 5 + M_corr] = z_alpha_corr;
  if (M_vcov > 0)
    regression_prior_raw[prior_start_alpha + 6 + M_corr : prior_start_alpha + 5 + M_corr + M_vcov] = z_alpha_vcov;
  if (n_marker_weight_means > 0)
    regression_prior_raw[prior_start_alpha + 6 + M_corr + M_vcov : prior_start_alpha + n_alpha_prior - 1] = z_marker_weight_mean;

  // Functional affine shifts are packed in exactly the same channel order as
  // their role vector in R: intercept then slope within each channel.
  if (n_iota_prior > 0) {
    int iota_position = prior_start_iota; // next free position in the common prior vector
    if (estimate_iota_intercept_cv > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cv - 1] = z_iota_intercept_cv; iota_position += estimate_iota_intercept_cv; }
    if (estimate_iota_slope_cv > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cv - 1] = z_iota_slope_cv; iota_position += estimate_iota_slope_cv; }
    if (estimate_iota_intercept_cs > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cs - 1] = z_iota_intercept_cs; iota_position += estimate_iota_intercept_cs; }
    if (estimate_iota_slope_cs > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cs - 1] = z_iota_slope_cs; iota_position += estimate_iota_slope_cs; }
    if (M_corr * estimate_iota_intercept_corr > 0) { regression_prior_raw[iota_position : iota_position + M_corr * estimate_iota_intercept_corr - 1] = z_iota_intercept_corr; iota_position += M_corr * estimate_iota_intercept_corr; }
    if (M_corr * estimate_iota_slope_corr > 0) { regression_prior_raw[iota_position : iota_position + M_corr * estimate_iota_slope_corr - 1] = z_iota_slope_corr; iota_position += M_corr * estimate_iota_slope_corr; }
    if (M_vcov * estimate_iota_intercept_vcov > 0) { regression_prior_raw[iota_position : iota_position + M_vcov * estimate_iota_intercept_vcov - 1] = z_iota_intercept_vcov; iota_position += M_vcov * estimate_iota_intercept_vcov; }
    if (M_vcov * estimate_iota_slope_vcov > 0) { regression_prior_raw[iota_position : iota_position + M_vcov * estimate_iota_slope_vcov - 1] = z_iota_slope_vcov; iota_position += M_vcov * estimate_iota_slope_vcov; }
    if (estimate_iota_intercept_cv_mean > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cv_mean - 1] = z_iota_intercept_cv_mean; iota_position += estimate_iota_intercept_cv_mean; }
    if (estimate_iota_slope_cv_mean > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cv_mean - 1] = z_iota_slope_cv_mean; iota_position += estimate_iota_slope_cv_mean; }
    if (estimate_iota_intercept_cv_marker > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cv_marker - 1] = z_iota_intercept_cv_marker; iota_position += estimate_iota_intercept_cv_marker; }
    if (estimate_iota_slope_cv_marker > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cv_marker - 1] = z_iota_slope_cv_marker; iota_position += estimate_iota_slope_cv_marker; }
    if (estimate_iota_intercept_cs_mean > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cs_mean - 1] = z_iota_intercept_cs_mean; iota_position += estimate_iota_intercept_cs_mean; }
    if (estimate_iota_slope_cs_mean > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cs_mean - 1] = z_iota_slope_cs_mean; iota_position += estimate_iota_slope_cs_mean; }
    if (estimate_iota_intercept_cs_marker > 0) { regression_prior_raw[iota_position : iota_position + estimate_iota_intercept_cs_marker - 1] = z_iota_intercept_cs_marker; iota_position += estimate_iota_intercept_cs_marker; }
    if (estimate_iota_slope_cs_marker > 0) regression_prior_raw[iota_position : iota_position + estimate_iota_slope_cs_marker - 1] = z_iota_slope_cs_marker;
  }

  if (P_vcov_sd > 0)
    regression_prior_raw[prior_start_vcov_sd : prior_start_vcov_sd + P_vcov_sd - 1] = vcov_sd_coefficient_raw;
  if (P_vcov_corr > 0)
    regression_prior_raw[prior_start_vcov_corr : prior_start_vcov_corr + P_vcov_corr - 1] = vcov_corr_coefficient_raw;
  if (K_event * p_w > 0) {
    for (event_cause in 1 : K_event) {
      int first_event_slope = prior_start_survival + (event_cause - 1) * p_w; // first covariate slope for this cause
      regression_prior_raw[first_event_slope : first_event_slope + p_w - 1] = z_gamma_w[event_cause];
    }
  }
  if (P_sigma > 0) regression_prior_raw[prior_start_sigma : prior_start_sigma + P_sigma - 1] = z_beta_sigma;
  if (P_nu > 0) regression_prior_raw[prior_start_nu : prior_start_nu + P_nu - 1] = z_beta_nu;
  if (P_phi > 0) regression_prior_raw[prior_start_phi : prior_start_phi + P_phi - 1] = z_beta_phi;
  if (P_alpha > 0) regression_prior_raw[prior_start_distributional_alpha : prior_start_distributional_alpha + P_alpha - 1] = z_beta_alpha;
  if (P_kappa > 0) regression_prior_raw[prior_start_kappa : prior_start_kappa + P_kappa - 1] = z_beta_kappa;
  if (P_tau > 0) regression_prior_raw[prior_start_tau : prior_start_tau + P_tau - 1] = z_beta_tau;

  regression_prior_effect = joinme_mixed_prior_transform(
    regression_prior_raw,
    prior_regression_family,
    prior_regression_mu,
    prior_regression_scale,
    prior_regression_horseshoe_local_index,
    prior_regression_horseshoe_group_index,
    horseshoe_local_regression,
    horseshoe_global_regression,
    horseshoe_slab_regression,
    prior_regression_horseshoe_slab_scale
  ); // all ordinary coefficients after their independently declared role priors

  if (P > 0) beta = regression_prior_effect[prior_start_beta : prior_start_beta + P - 1];
  alpha_prior_effect = regression_prior_effect[prior_start_alpha : prior_start_alpha + n_alpha_prior - 1];
  if (K_event * p_w > 0)
    for (event_cause in 1 : K_event) {
      int first_event_slope = prior_start_survival + (event_cause - 1) * p_w; // first transformed covariate slope for this cause
      gamma_w[event_cause] = regression_prior_effect[first_event_slope : first_event_slope + p_w - 1];
    }
  if (P_sigma > 0) beta_sigma = regression_prior_effect[prior_start_sigma : prior_start_sigma + P_sigma - 1];
  if (P_nu > 0) beta_nu = regression_prior_effect[prior_start_nu : prior_start_nu + P_nu - 1];
  if (P_phi > 0) beta_phi = regression_prior_effect[prior_start_phi : prior_start_phi + P_phi - 1];
  if (P_alpha > 0) beta_alpha = regression_prior_effect[prior_start_distributional_alpha : prior_start_distributional_alpha + P_alpha - 1];
  if (P_kappa > 0) beta_kappa = regression_prior_effect[prior_start_kappa : prior_start_kappa + P_kappa - 1];
  if (P_tau > 0) beta_tau = regression_prior_effect[prior_start_tau : prior_start_tau + P_tau - 1];

  vector[M_cov] alpha_L = rep_vector(0, M_cov); // packed SD and correlation intercepts retained in the established lower-triangular reporting order
  array[Q_idm] vector[K_cov_sd] beta_L_sd; // slopes unique to each standard-deviation coordinate and its formulaVCov$sd design
  array[M_corr] vector[K_cov_corr] beta_L_corr; // slopes unique to each off-diagonal partial correlation and its formulaVCov$corr design

  // STEP 2: restore the two covariance-regression coefficient blocks.
  //
  // Each block has its own formula, dimension and prior family. Packing all
  // intercepts before the row-major slope matrix gives R, simulation and Stan
  // one deterministic parameter order even when the two formulae differ.
  {
    vector[P_vcov_sd] vcov_sd_coefficient = rep_vector(0, P_vcov_sd); // transformed SD intercepts and slopes
    vector[P_vcov_corr] vcov_corr_coefficient = rep_vector(0, P_vcov_corr); // transformed partial-correlation intercepts and slopes
    int correlation_coordinate = 1; // next row-major off-diagonal coordinate when reconstructing packed alpha_L

    if (P_vcov_sd > 0)
      vcov_sd_coefficient = regression_prior_effect[prior_start_vcov_sd : prior_start_vcov_sd + P_vcov_sd - 1];
    if (P_vcov_corr > 0)
      vcov_corr_coefficient = regression_prior_effect[prior_start_vcov_corr : prior_start_vcov_corr + P_vcov_corr - 1];

    if (K_cov_sd > 0) {
      for (r in 1 : Q_idm) {
        int first_slope = Q_idm + (r - 1) * K_cov_sd + 1; // first packed SD slope for row r
        int final_slope = Q_idm + r * K_cov_sd; // final packed SD slope for row r
        beta_L_sd[r] = vcov_sd_coefficient[first_slope : final_slope];
      }
    }
    if (K_cov_corr > 0) {
      for (correlation in 1 : M_corr) {
        int first_slope = M_corr + (correlation - 1) * K_cov_corr + 1; // first packed slope for this off-diagonal coordinate
        int final_slope = M_corr + correlation * K_cov_corr; // final packed slope for this off-diagonal coordinate
        beta_L_corr[correlation] = vcov_corr_coefficient[first_slope : final_slope];
      }
    }
    for (m in 1 : M_cov) {
      if (r_idx[m] == c_idx[m]) {
        alpha_L[m] = vcov_sd_coefficient[r_idx[m]];
      } else {
        alpha_L[m] = vcov_corr_coefficient[correlation_coordinate];
        correlation_coordinate += 1;
      }
    }
  }

