  /**
   * @brief Time-scaled coefficient vectors and derived random effects for both longitudinal and survival models.
   *
   * This block constructs:
   * 1. Scaled versions of fixed and random effect SDs where time transformations are applied
   * 2. Cholesky factors for random effect covariance structures
   * 3. Marker-average effects used in survival model
   * 4. ID-specific variance structure for covariance regression
   */

  // -------------------- Effective (scaled-time) coefficient vectors for use with scaled design matrices
  vector[P] beta_scaled = beta;
  if (n_time_beta > 0) {
    for (k in 1 : n_time_beta) 
      beta_scaled[idx_time_beta[k]] = beta[idx_time_beta[k]] * tmax;
  }
  
  // Build scale vectors for SDs so random slope components are not inadvertently shrunk
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
  
  vector[Q_idm] tau_w_scaled = tau_w;
  if (n_time_widm > 0) {
    for (k in 1 : n_time_widm) 
      tau_w_scaled[idx_time_widm[k]] = tau_w[idx_time_widm[k]] * tmax;
  }
  
  // -------------------- id RE covariance
  matrix[R_id, R_id] L_u;
  if (indep_id_re == 1) 
    L_u = diag_matrix(tau_u_scaled);
  else 
    L_u = diag_pre_multiply(tau_u_scaled, Lcorr_u);
  
  array[n_id] vector[R_id] u_id;
  for (i in 1 : n_id) 
    u_id[i] = L_u * z_u[i];
  
  // -------------------- marker-only covariance (only meaningful if R_mk>0)
  matrix[R_mk, R_mk] L_v;
  array[D] vector[R_mk] v_marker;
  if (R_mk > 0) {
    if (indep_marker_re == 1) 
      L_v = diag_matrix(tau_v_scaled);
    else 
      L_v = diag_pre_multiply(tau_v_scaled, Lcorr_v);
    
    for (d in 1 : D) 
      v_marker[d] = L_v * z_v[d];
  }
  
  // -------------------- marker-by-id latent covariance
  matrix[Q_idm, Q_idm] L_w;
  if (indep_marker_byid_latent_re == 1) 
    L_w = diag_matrix(tau_w_scaled);
  else 
    L_w = diag_pre_multiply(tau_w_scaled, Lcorr_w);
  
  // Construct z_w[i,d], optionally cross-correlated with v_marker[d]
  array[n_id, D] vector[Q_idm] z_w;
  for (i in 1 : n_id) {
    for (d in 1 : D) {
      vector[Q_idm] cross = rep_vector(0.0, Q_idm);
      if (allow_marker_crosscorr == 1 && R_mk > 0)
        cross = B_cross * v_marker[d];
      z_w[i, d] = cross + L_w * z_w_lat[i, d];
    }
  }
  
  // -------------------- id random effect for covariance regression
  vector[n_id] u_L;
  for (i in 1 : n_id) 
    u_L[i] = tau_L * z_L[i];
  
  // -------------------- id-specific Cholesky factors L_i
  array[n_id] matrix[Q_idm, Q_idm] L_i;
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] Li = rep_matrix(0.0, Q_idm, Q_idm);
    for (m in 1 : M_cov) {
      int r = r_idx[m];
      int c = c_idx[m];
      real lp = alpha_L[m] + dot_product(beta_L[m], to_vector(Xcov[i]'))
                + lambda_L[m] * u_L[i];
      Li[r, c] = (r == c) ? log1p_exp(lp) : lp;
    }
    L_i[i] = Li;
  }
  
  // -------------------- marker-by-id scaled effects: w_idscaled[i,d] = L_i[i] * z_w[i,d]
  array[n_id, D] vector[Q_idm] w_idscaled;
  for (i in 1 : n_id)
    for (d in 1 : D)
      w_idscaled[i, d] = L_i[i] * z_w[i, d];
  
  // -------------------- marker weights (optional shrinkage around base weights)
  // If estimate_marker_weights = 1, use a logistic-normal perturbation of the
  // base weights with scale tau_marker_weights. Otherwise use marker_weights as-is.
  vector[D] marker_weights_eff;
  {
    vector[D] base_log = log(marker_weights + 1e-12);
    if (estimate_marker_weights == 1) {
      marker_weights_eff = softmax(base_log + tau_marker_weights * z_marker_weights);
    } else {
      marker_weights_eff = marker_weights;
    }
  }

  // -------------------- marker averages for survival association (weighted mean across markers)
  // These averages are used for both current value (CV) and current slope (CS) terms.
  real sum_w = sum(marker_weights_eff);
  if (sum_w <= 0) sum_w = 1;
  vector[R_mk] vbar;
  if (R_mk > 0) {
    for (r in 1 : R_mk) {
      real acc = 0;
      for (d in 1 : D) 
        acc += marker_weights_eff[d] * v_marker[d][r];
      vbar[r] = acc / sum_w;
    }
  }
  
  array[n_id] vector[Q_idm] zbar;
  for (i in 1 : n_id) {
    for (q in 1 : Q_idm) {
      real acc = 0;
      for (d in 1 : D)
        acc += marker_weights_eff[d] * z_w[i, d][q];
      zbar[i][q] = acc / sum_w;
    }
  }
  
  array[n_id] vector[Q_idm] wbar_i;
  for (i in 1 : n_id)
    wbar_i[i] = L_i[i] * zbar[i];
  
  // -------------------- Effective association coefficients (flags applied)
  real a_cv_total = assoc_cv_total * alpha_cv_total;
  real a_cs_total = assoc_cs_total * alpha_cs_total;
  real a_cv_mean = assoc_cv_mean * alpha_cv_mean;
  real a_cv_marker = assoc_cv_marker * alpha_cv_marker;
  real a_cs_mean = assoc_cs_mean * alpha_cs_mean;
  real a_cs_marker = assoc_cs_marker * alpha_cs_marker;
  vector[Q_idm] a_vcov_var = assoc_vcov * alpha_vcov_var;
