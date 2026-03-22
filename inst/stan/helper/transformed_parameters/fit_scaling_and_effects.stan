  /**
  * @brief Time-scaled coefficient vectors and derived random effects for both longitudinal and survival models.
   *
   * This block constructs:
   * 1. Scaled versions of fixed and random effect SDs where time transformations are applied
   * 2. Cholesky factors for random effect covariance structures
   * 3. Marker-average effects used in survival model
   * 4. ID-specific variance structure for covariance regression
   */

  /* -------------------- Effective (scaled-time) coefficient vectors for use with scaled design matrices */
  vector[P] beta_scaled = beta;
  if (n_time_beta > 0) {
    for (k in 1 : n_time_beta) 
      beta_scaled[idx_time_beta[k]] = beta[idx_time_beta[k]] * tmax;
  }
  
  /* Build scale vectors for SDs so random slope components are not inadvertently shrunk */
  vector[R_id] tau_u_scaled = tau_u;
  if (n_time_uid > 0) {
    for (k in 1 : n_time_uid) 
      tau_u_scaled[idx_time_uid[k]] = tau_u[idx_time_uid[k]] * tmax;
  }
  
  vector[R_mk] tau_v_scaled;
  if (R_mk > 0) {
    tau_v_scaled = tau_v;
    if (n_time_vmk > 0) {
      // idx_time_vmk are validated implicitly by use; if out-of-range it will error at runtime.
      for (k in 1 : n_time_vmk) 
        tau_v_scaled[idx_time_vmk[k]] = tau_v[idx_time_vmk[k]] * tmax;
    }
  }

  // Marker-by-id covariance regression is parameterised on the original-time
  // scale. We only row-scale time-indexed marker-by-id coefficients when they
  // are paired with the scaled-time longitudinal design matrices.
  vector[Q_idm] row_scale_idm = rep_vector(1.0, Q_idm);
  if (n_time_idm > 0) {
    for (k in 1 : n_time_idm)
      row_scale_idm[idx_time_idm[k]] = tmax;
  }
  
  /* -------------------- id RE covariance */
  matrix[R_id, R_id] L_u; // Cholesky factor for id REs on scaled time
  if (indep_id_re == 1) 
    L_u = diag_matrix(tau_u_scaled);
  else 
    L_u = diag_pre_multiply(tau_u_scaled, Lcorr_u);
  
  array[n_id] vector[R_id] u_id; // realized id random effects
  for (i in 1 : n_id) 
    u_id[i] = L_u * z_u[i];
  
  /* -------------------- marker-only covariance (only meaningful if R_mk>0) */
  matrix[R_mk, R_mk] L_v; // Cholesky factor for marker REs
  array[D] vector[R_mk] v_marker; // realized marker random effects
  if (R_mk > 0) {
    if (indep_marker_re == 1) 
      L_v = diag_matrix(tau_v_scaled);
    else 
      L_v = diag_pre_multiply(tau_v_scaled, Lcorr_v);
    
    for (d in 1 : D) 
      v_marker[d] = L_v * z_v[d];
  }
  
  /* -------------------- marker-by-id latent effects */
  // The marker-by-id latent seeds are now always iid standard normal.
  // Baseline covariance and subject-specific covariance changes are both
  // encoded by L_i below, which avoids a redundant global covariance layer.
  // Optional cross-correlation with marker-level effects is still added here
  // before the subject-specific covariance scaling is applied.
  array[n_id, D] vector[Q_idm] z_w; // latent marker-id effects after optional cross-correlation
  for (i in 1 : n_id) {
    for (d in 1 : D) {
      vector[Q_idm] cross = rep_vector(0.0, Q_idm);
      if (allow_marker_crosscorr == 1 && R_mk > 0)
        cross = B_cross * v_marker[d];
      z_w[i, d] = cross + z_w_lat[i, d];
    }
  }
  
  /* -------------------- id-specific Cholesky factors L_i */
  array[n_id] matrix[Q_idm, Q_idm] L_i; // subject-specific Cholesky factor
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] Li = rep_matrix(0.0, Q_idm, Q_idm); // local accumulator
    for (m in 1 : M_cov) {
      int r = r_idx[m];
      int c = c_idx[m];
      real lp = alpha_L[m] + dot_product(beta_L[m], to_vector(Xcov[i]'))
        + lambda_L[m] * z_L[i][m]; // linear predictor for L_i element
      if (r == c) {
        Li[r, c] = (vcov_diag_link == 1) ? exp(lp) : log1p_exp(lp);
      } else {
        Li[r, c] = lp;
      }
    }
    L_i[i] = Li;
  }
  
  /* -------------------- marker-by-id scaled effects: w_idm[i,d] = L_i_eff[i] * z_w[i,d] */
  array[n_id, D] vector[Q_idm] w_idm; // scaled marker-id effects
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] L_i_eff = diag_pre_multiply(row_scale_idm, L_i[i]);
    for (d in 1 : D)
      w_idm[i, d] = L_i_eff * z_w[i, d];
  }
  
  /* -------------------- marker weights (signed) */
  // Goal:
  // - Allow positive and negative marker contributions.
  // - Keep perturbation model simple and directly interpretable.
  // - No normalisation is applied; raw weights are used directly.
  // Rule:
  // - Start from signed base weights (provided from standata) as a prior offset.
  // - If estimation is enabled, add a signed standard-normal perturbation z.
  // - Map to effective weights via w_raw (no normalisation).
  // - Interpret marker_weights_eff as marker-intensity multipliers for association.
  vector[D] marker_weights_eff; // effective signed marker weights
  {
    if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1) {
      marker_weights_eff = marker_weights + z_marker_weights;
    } else {
      marker_weights_eff = marker_weights;
    }
  }

  /* -------------------- effective spline coefficients for transforms */
  // For Stan-estimated penalised splines, we estimate only the SHAPE here.
  // We anchor the transform so:
  // - the first coefficient is 0,
  // - the last coefficient is 1,
  // - all intermediate coefficients are monotone increasing.
  // This avoids confounding transform scale with the association coefficient.
  // Implementation detail:
  // - we build positive increments with `softmax()`,
  // - those increments sum to 1 exactly,
  // - cumulative sums therefore give a stable unit-span monotone curve
  //   without any divide-by-a-nearly-zero normalisation step.
  vector[n_coeff_cv] coeff_cv_eff = coeff_cv;
  vector[n_coeff_cs] coeff_cs_eff = coeff_cs;
  matrix[M_corr, n_coeff_corr] coeff_corr_eff = coeff_corr;
  matrix[M_vcov, n_coeff_vcov] coeff_vcov_eff = coeff_vcov;
  vector[n_coeff_cv_mean] coeff_cv_mean_eff = coeff_cv_mean;
  vector[n_coeff_cv_marker] coeff_cv_marker_eff = coeff_cv_marker;
  vector[n_coeff_cs_mean] coeff_cs_mean_eff = coeff_cs_mean;
  vector[n_coeff_cs_marker] coeff_cs_marker_eff = coeff_cs_marker;

  if (estimate_spline_cv == 1 && n_coeff_cv > 1 && n_free_spline_cv == n_coeff_cv - 1) {
    vector[n_coeff_cv - 1] delta = softmax(z_spline_cv);
    coeff_cv_eff[1] = 0;
    for (j in 2:n_coeff_cv) coeff_cv_eff[j] = coeff_cv_eff[j - 1] + delta[j - 1];
  }
  if (estimate_spline_cs == 1 && n_coeff_cs > 1 && n_free_spline_cs == n_coeff_cs - 1) {
    vector[n_coeff_cs - 1] delta = softmax(z_spline_cs);
    coeff_cs_eff[1] = 0;
    for (j in 2:n_coeff_cs) coeff_cs_eff[j] = coeff_cs_eff[j - 1] + delta[j - 1];
  }
  if (estimate_spline_corr == 1 && n_coeff_corr > 1 && n_free_spline_corr == n_coeff_corr - 1) {
    for (m in 1:M_corr) {
      vector[n_coeff_corr - 1] delta = softmax(to_vector(row(z_spline_corr, m)));
      coeff_corr_eff[m, 1] = 0;
      for (j in 2:n_coeff_corr) coeff_corr_eff[m, j] = coeff_corr_eff[m, j - 1] + delta[j - 1];
    }
  }
  if (estimate_spline_vcov == 1 && n_coeff_vcov > 1 && n_free_spline_vcov == n_coeff_vcov - 1) {
    for (m in 1:M_vcov) {
      vector[n_coeff_vcov - 1] delta = softmax(to_vector(row(z_spline_vcov, m)));
      coeff_vcov_eff[m, 1] = 0;
      for (j in 2:n_coeff_vcov) coeff_vcov_eff[m, j] = coeff_vcov_eff[m, j - 1] + delta[j - 1];
    }
  }
  if (estimate_spline_cv_mean == 1 && n_coeff_cv_mean > 1 && n_free_spline_cv_mean == n_coeff_cv_mean - 1) {
    vector[n_coeff_cv_mean - 1] delta = softmax(z_spline_cv_mean);
    coeff_cv_mean_eff[1] = 0;
    for (j in 2:n_coeff_cv_mean) coeff_cv_mean_eff[j] = coeff_cv_mean_eff[j - 1] + delta[j - 1];
  }
  if (estimate_spline_cv_marker == 1 && n_coeff_cv_marker > 1 && n_free_spline_cv_marker == n_coeff_cv_marker - 1) {
    vector[n_coeff_cv_marker - 1] delta = softmax(z_spline_cv_marker);
    coeff_cv_marker_eff[1] = 0;
    for (j in 2:n_coeff_cv_marker) coeff_cv_marker_eff[j] = coeff_cv_marker_eff[j - 1] + delta[j - 1];
  }
  if (estimate_spline_cs_mean == 1 && n_coeff_cs_mean > 1 && n_free_spline_cs_mean == n_coeff_cs_mean - 1) {
    vector[n_coeff_cs_mean - 1] delta = softmax(z_spline_cs_mean);
    coeff_cs_mean_eff[1] = 0;
    for (j in 2:n_coeff_cs_mean) coeff_cs_mean_eff[j] = coeff_cs_mean_eff[j - 1] + delta[j - 1];
  }
  if (estimate_spline_cs_marker == 1 && n_coeff_cs_marker > 1 && n_free_spline_cs_marker == n_coeff_cs_marker - 1) {
    vector[n_coeff_cs_marker - 1] delta = softmax(z_spline_cs_marker);
    coeff_cs_marker_eff[1] = 0;
    for (j in 2:n_coeff_cs_marker) coeff_cs_marker_eff[j] = coeff_cs_marker_eff[j - 1] + delta[j - 1];
  }

  /* -------------------- marker averages for survival association (unweighted mean across markers) */
  // For transformed marker-aggregated CV terms, signed marker weights are applied
  // only after transformation. vbar/wbar_i stay unweighted means to avoid
  // double-weighting when transformed marker aggregation is enabled.
  vector[R_mk] vbar; // unweighted mean of marker REs
  if (R_mk > 0) {
    for (r in 1 : R_mk) {
      real acc = 0;
      for (d in 1 : D) 
        acc += v_marker[d][r];
      vbar[r] = acc / D;
    }
  }
  
  array[n_id] vector[Q_idm] zbar; // unweighted mean of marker-id latents
  for (i in 1 : n_id) {
    for (q in 1 : Q_idm) {
      real acc = 0;
      for (d in 1 : D)
        acc += z_w[i, d][q];
      zbar[i][q] = acc / D;
    }
  }
  
  array[n_id] vector[Q_idm] wbar_i; // scaled weighted mean per subject
  for (i in 1 : n_id)
    wbar_i[i] = diag_pre_multiply(row_scale_idm, L_i[i]) * zbar[i];
  
  /* -------------------- Effective association coefficients (flags applied) */
  // Non-centred association construction with flexible sign constraints:
  // if any marker-weight-involving association is active, only the first
  // active latent in the order (cv_total, cs_total, cv_marker, cs_marker)
  // is constrained positive; all others remain unconstrained.
  real z_alpha_cv_total_eff = z_alpha_cv_total;
  real z_alpha_cs_total_eff = z_alpha_cs_total;
  real z_alpha_cv_marker_eff = z_alpha_cv_marker;
  real z_alpha_cs_marker_eff = z_alpha_cs_marker;

  if ((assoc_cv_total + assoc_cs_total + assoc_cv_marker + assoc_cs_marker) > 0) {
    if (assoc_cv_total == 1) {
      z_alpha_cv_total_eff = abs(z_alpha_cv_total);
    } else if (assoc_cs_total == 1) {
      z_alpha_cs_total_eff = abs(z_alpha_cs_total);
    } else if (assoc_cv_marker == 1) {
      z_alpha_cv_marker_eff = abs(z_alpha_cv_marker);
    } else if (assoc_cs_marker == 1) {
      z_alpha_cs_marker_eff = abs(z_alpha_cs_marker);
    }
  }

  real alpha_cv_total = z_alpha_cv_total_eff * sd_alpha_cv_total;
  real alpha_cs_total = z_alpha_cs_total_eff * sd_alpha_cs_total;
  real alpha_cv_mean = z_alpha_cv_mean * sd_alpha_cv_mean;
  real alpha_cs_mean = z_alpha_cs_mean * sd_alpha_cs_mean;
  real alpha_cv_marker = z_alpha_cv_marker_eff * sd_alpha_cv_marker;
  real alpha_cs_marker = z_alpha_cs_marker_eff * sd_alpha_cs_marker;
  vector[M_corr] alpha_corr_eff = s_corr * alpha_corr; // effective corr coefficients on the model scale
  vector[M_vcov] alpha_vcov_eff = s_vcov * alpha_vcov; // effective vcov coefficients on the model scale

  real a_cv_total = assoc_cv_total * alpha_cv_total; // total CV coefficient
  real a_cs_total = assoc_cs_total * alpha_cs_total; // total CS coefficient
  real a_cv_mean = assoc_cv_mean * alpha_cv_mean;    // mean CV coefficient
  real a_cv_marker = assoc_cv_marker * alpha_cv_marker; // marker CV coefficient
  real a_cs_mean = assoc_cs_mean * alpha_cs_mean;    // mean CS coefficient
  real a_cs_marker = assoc_cs_marker * alpha_cs_marker; // marker CS coefficient
  vector[M_corr] a_corr = assoc_corr * alpha_corr_eff; // corr correlation coefficients actually used in the hazard
  vector[M_vcov] a_vcov = assoc_vcov * alpha_vcov_eff; // vcov association coefficients actually used in the hazard
