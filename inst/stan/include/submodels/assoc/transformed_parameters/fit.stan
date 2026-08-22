/**
 * @file assoc/transformed_parameters/fit.stan
 * @brief Construct association summaries and active hazard coefficients.
 *
 * @details
 * Marker-level random effects are averaged for mean channels, sign conventions are applied where estimated weights require identification, and inactive association channels are set to zero before the event likelihood is evaluated.
 */
  /* -------------------- marker averages for survival association (unweighted mean across markers) */
  // For transformed marker-aggregated CV terms, signed marker weights are applied
  // only after transformation. vbar/wbar_i stay unweighted means to avoid
  // double-weighting when transformed marker aggregation is enabled.
  vector[R_mk] vbar; // unweighted mean of marker REs
  if (R_mk > 0) {
    for (r in 1 : R_mk) {
      real acc = 0; // Role: running accumulator.
      for (d in 1 : D) 
        acc += v_marker[d][r];
      vbar[r] = acc / D;
    }
  }
  
  array[n_id] vector[Q_idm] zbar; // unweighted mean of marker-id latents
  for (i in 1 : n_id) {
    for (q in 1 : Q_idm) {
      real acc = 0; // Role: running accumulator.
      for (d in 1 : D)
        acc += z_w[i, d][q];
      zbar[i][q] = acc / D;
    }
  }
  
  array[n_id] vector[Q_idm] wbar_i; // scaled weighted mean per subject
  for (i in 1 : n_id)
    wbar_i[i] = L_i[i] * zbar[i];
  
  /* -------------------- Effective association coefficients (flags applied) */
  // Weighted marker features have a sign symmetry: reversing every marker
  // weight and the corresponding association coefficient leaves the event
  // predictor unchanged.  Constraining the association coefficient to be
  // non-negative chooses one of these equivalent orientations without
  // changing its magnitude.  The common marker-weight location remains part
  // of every effective weight and is not replaced by this convention.
  real z_alpha_cv_total_constrained = alpha_prior_effect[1]; // transformed latent for total current-value association
  real z_alpha_cs_total_constrained = alpha_prior_effect[2]; // transformed latent for total current-slope association
  real z_alpha_cv_marker_constrained = alpha_prior_effect[5]; // transformed latent for marker current-value association
  real z_alpha_cs_marker_constrained = alpha_prior_effect[6]; // transformed latent for marker current-slope association

  if (assoc_cv_total == 1) z_alpha_cv_total_constrained = abs(alpha_prior_effect[1]);
  if (assoc_cs_total == 1) z_alpha_cs_total_constrained = abs(alpha_prior_effect[2]);
  if (assoc_cv_marker == 1) z_alpha_cv_marker_constrained = abs(alpha_prior_effect[5]);
  if (assoc_cs_marker == 1) z_alpha_cs_marker_constrained = abs(alpha_prior_effect[6]);

  /*
   * These quantities are the association coefficients used by the hazard.
   * They are not multiplied by an additional unidentified random scale: the
   * location and scale supplied in `jm_priors(assoc = list(slope = ...))` therefore describe
   * the coefficient itself. The suffix "scaled" is retained because the
   * downstream likelihood has historically used these names.
   */
  real alpha_cv_total_scaled = z_alpha_cv_total_constrained; // Role: total current-value association coefficient on the hazard scale.
  real alpha_cs_total_scaled = z_alpha_cs_total_constrained; // Role: total current-slope association coefficient on the hazard scale.
  real alpha_cv_mean_scaled = alpha_prior_effect[3]; // Role: subject-mean current-value association coefficient on the hazard scale.
  real alpha_cs_mean_scaled = alpha_prior_effect[4]; // Role: subject-mean current-slope association coefficient on the hazard scale.
  real alpha_cv_marker = z_alpha_cv_marker_constrained; // Role: marker current-value association coefficient on the hazard scale.
  real alpha_cs_marker = z_alpha_cs_marker_constrained; // Role: marker current-slope association coefficient on the hazard scale.
  vector[M_corr] alpha_corr_scaled; // transformed correlation-association coefficients
  vector[M_vcov] alpha_vcov_scaled; // transformed covariance-association coefficients
  if (M_corr > 0) alpha_corr_scaled = alpha_prior_effect[7 : 6 + M_corr];
  if (M_vcov > 0) alpha_vcov_scaled = alpha_prior_effect[7 + M_corr : 6 + M_corr + M_vcov];

  real a_cv_total = assoc_cv_total * alpha_cv_total_scaled; // total CV coefficient
  real a_cs_total = assoc_cs_total * alpha_cs_total_scaled; // total CS coefficient
  real a_cv_mean = assoc_cv_mean * alpha_cv_mean_scaled;    // mean CV coefficient
  real a_cv_marker = assoc_cv_marker * alpha_cv_marker; // marker CV coefficient
  real a_cs_mean = assoc_cs_mean * alpha_cs_mean_scaled;    // mean CS coefficient
  real a_cs_marker = assoc_cs_marker * alpha_cs_marker; // marker CS coefficient
  vector[M_corr] a_corr = assoc_corr * alpha_corr_scaled; // corr coefficients used in the hazard
  vector[M_vcov] a_vcov = assoc_vcov * alpha_vcov_scaled; // vcov coefficients used in the hazard
