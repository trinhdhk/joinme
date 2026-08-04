  /**
  * @brief Time-scaled coefficient vectors and derived random effects for both longitudinal and survival models.
   *
   * This block constructs:
   * 1. Scaled versions of fixed and random effect SDs where time transformations are applied
   * 2. Cholesky factors for random effect covariance structures
   * 3. Marker-average effects used in survival model
  * 4. ID-specific variance structure for covariance regression
  */

  vector[M_cov] alpha_L = rep_vector(0, M_cov); // packed SD and correlation intercepts retained in the established lower-triangular reporting order
  array[Q_idm] vector[K_cov_sd] beta_L_sd; // slopes unique to each standard-deviation coordinate and its formulaVCov$sd design
  array[M_corr] vector[K_cov_corr] beta_L_corr; // slopes unique to each off-diagonal partial correlation and its formulaVCov$corr design

  // STEP 1: transform the two covariance-regression coefficient blocks.
  //
  // Each block has its own formula, dimension and prior family. Packing all
  // intercepts before the row-major slope matrix gives R, simulation and Stan
  // one deterministic parameter order even when the two formulae differ.
  {
    vector[P_vcov_sd] vcov_sd_coefficient = joinme_prior_transform(
      vcov_sd_coefficient_raw,
      prior_vcov_sd_family,
      prior_vcov_sd_mu,
      prior_vcov_sd_scale,
      horseshoe_local_vcov_sd,
      horseshoe_global_vcov_sd,
      horseshoe_slab_vcov_sd,
      prior_vcov_sd_slab_scale
    ); // SD intercepts and slopes after their selected prior transformation
    vector[P_vcov_corr] vcov_corr_coefficient = joinme_prior_transform(
      vcov_corr_coefficient_raw,
      prior_vcov_corr_family,
      prior_vcov_corr_mu,
      prior_vcov_corr_scale,
      horseshoe_local_vcov_corr,
      horseshoe_global_vcov_corr,
      horseshoe_slab_vcov_corr,
      prior_vcov_corr_slab_scale
    ); // off-diagonal correlation intercepts and slopes after their independent transformation
    int correlation_coordinate = 1; // next row-major off-diagonal coordinate when reconstructing packed alpha_L

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

  /* -------------------- Effective (scaled-time) coefficient vectors for use with scaled design matrices */
  vector[P] beta_scaled = beta; // Role: coefficient scaled.
  if (n_time_beta > 0) {
    for (k in 1 : n_time_beta) 
      beta_scaled[idx_time_beta[k]] = beta[idx_time_beta[k]] * tmax;
  }
  
  /* Build scale vectors for SDs so random slope components are not inadvertently shrunk */
  vector[R_id] tau_u_scaled = tau_u; // Role: quantile u scaled.
  if (n_time_uid > 0) {
    for (k in 1 : n_time_uid) 
      tau_u_scaled[idx_time_uid[k]] = tau_u[idx_time_uid[k]] * tmax;
  }
  
  vector[R_mk] tau_v_scaled; // Role: quantile v scaled.
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
  vector[Q_idm] row_scale_idm = rep_vector(1.0, Q_idm); // Role: row scale idm.
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
    vector[D * R_mk] marker_raw; // packed standardised marker parameters in marker-major order
    vector[D * R_mk] marker_prior_effect; // packed marker parameters after their selected unit-scale prior transformation
    for (d in 1 : D) {
      int marker_start = (d - 1) * R_mk + 1; // first packed coordinate for marker d
      int marker_end = d * R_mk; // final packed coordinate for marker d
      marker_raw[marker_start : marker_end] = z_v[d];
    }
    marker_prior_effect = joinme_prior_transform(
      marker_raw,
      prior_marker_family,
      prior_marker_mu,
      prior_marker_scale,
      horseshoe_local_marker,
      horseshoe_global_marker,
      horseshoe_slab_marker,
      prior_marker_slab_scale
    );
    if (indep_marker_re == 1)
      L_v = diag_matrix(tau_v_scaled);
    else 
      L_v = diag_pre_multiply(tau_v_scaled, Lcorr_v);
    
    for (d in 1 : D) {
      int marker_start = (d - 1) * R_mk + 1; // first transformed coordinate for marker d
      int marker_end = d * R_mk; // final transformed coordinate for marker d
      v_marker[d] = L_v * marker_prior_effect[marker_start : marker_end];
    }
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
      vector[Q_idm] cross = rep_vector(0.0, Q_idm); // Role: cross term.
      if (allow_marker_crosscorr == 1 && R_mk > 0)
        cross = B_cross * v_marker[d];
      z_w[i, d] = cross + z_w_lat[i, d];
    }
  }
  
  /* -------------------- id-specific Cholesky factors L_i = SD_i * K_i */
  array[n_id] matrix[Q_idm, Q_idm] L_i; // subject-specific covariance Cholesky factor
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] Li = rep_matrix(0.0, Q_idm, Q_idm); // local accumulator
    vector[Q_idm] sd_i = rep_vector(1.0, Q_idm); // subject-specific standard deviations
    int m_pos = 1; // Role: m pos.

    // Step 1: reconstruct the subject-specific SD regression on the diagonal.
    //
    // lambda_L[m] is a scalar residual-heterogeneity loading for exactly one
    // packed covariance coordinate m. Since z_L[i][m] has unit scale in the
    // ordinary model, lambda_L[m]^2 is the conditional variance contributed
    // to this linear predictor before applying exp or softplus. There are no
    // cross-coordinate loadings: in vector notation this term is
    // diag_matrix(lambda_L) * z_L[i], not a dense matrix times z_L[i].
    for (m in 1 : M_cov) {
      int r = r_idx[m]; // Role: row position.
      int c = c_idx[m]; // Role: column position.
      if (r == c) {
        real lp = alpha_L[m] + dot_product(beta_L_sd[r], to_vector(Xcov_sd[i]'))
          + lambda_L[m] * z_L[i][m];
        sd_i[r] = (vcov_diag_link == 1) ? exp(lp) : log1p_exp(lp);
      }
    }

    // Step 2: rebuild the Cholesky-correlation rows from tanh-linked
    // row-specific partial correlations. This keeps K_i valid for every subject.
    // For an off-diagonal coordinate, lambda_L[m] is the residual scale before
    // tanh; it is therefore not itself a correlation standard deviation.
    for (r in 1 : Q_idm) {
      if (r == 1) {
        Li[1, 1] = sd_i[1];
        if (indep_idmarker_cov == 0)
          m_pos += 1;
      } else if (indep_idmarker_cov == 1) {
        Li[r, r] = sd_i[r];
        m_pos += 1;
      } else {
        real scale_prod = 1.0; // Role: scale product.
        for (c in 1 : r) {
          if (c < r) {
            int correlation_coordinate = ((r - 1) * (r - 2)) %/% 2 + c; // row-major off-diagonal index: (2,1), (3,1), (3,2), ...
            real lp = alpha_L[m_pos] + dot_product(beta_L_corr[correlation_coordinate], to_vector(Xcov_corr[i]'))
              + lambda_L[m_pos] * z_L[i][m_pos];
            real z_rc = tanh(lp); // Role: standardised latent value rc.
            Li[r, c] = sd_i[r] * scale_prod * z_rc;
            scale_prod *= sqrt(fmax(1e-12, 1.0 - square(z_rc)));
          } else {
            Li[r, r] = sd_i[r] * scale_prod;
          }
          m_pos += 1;
        }
      }
    }

    L_i[i] = Li;
  }
  
  /* -------------------- marker-by-id scaled effects: w_idm[i,d] = L_i_eff[i] * z_w[i,d] */
  array[n_id, D] vector[Q_idm] w_idm; // scaled marker-id effects
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] L_i_eff = diag_pre_multiply(row_scale_idm, L_i[i]); // Role: Cholesky factor i eff.
    for (d in 1 : D)
      w_idm[i, d] = L_i_eff * z_w[i, d];
  }
  
  /* -------------------- marker weights (signed) */
  // Goal:
  // - Allow positive and negative marker contributions.
  // - Support either one shared marker-weight structure across all weighted
  //   marker-based association terms or one structure per active term.
  // - Keep perturbation model directly interpretable by adding signed latent
  //   perturbations to the supplied base weights.
  matrix[n_marker_weight_sets, D] z_marker_weight_sets = rep_matrix(0.0, n_marker_weight_sets, D); // Role: standardised latent value marker weight sets.
  vector[D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets] marker_weight_prior_effect = joinme_prior_transform(
    z_marker_weights,
    prior_marker_weight_family,
    prior_marker_weight_mu,
    prior_marker_weight_scale,
    horseshoe_local_marker_weight,
    horseshoe_global_marker_weight,
    horseshoe_slab_marker_weight,
    prior_marker_weight_slab_scale
  ); // marker-weight perturbations after the selected fixed-unit-scale family transformation
  vector[D] marker_weights_eff_cv_total = marker_weights_cv_total; // Role: marker weights eff current value total.
  vector[D] marker_weights_eff_cs_total = marker_weights_cs_total; // Role: marker weights eff current slope total.
  vector[D] marker_weights_eff_cv_marker = marker_weights_cv_marker; // Role: marker weights eff current value marker.
  vector[D] marker_weights_eff_cs_marker = marker_weights_cs_marker; // Role: marker weights eff current slope marker.
  {
    if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1) {
      for (s in 1:n_marker_weight_sets) {
        int start_pos = (s - 1) * D + 1; // Role: starting pos.
        int end_pos = s * D; // Role: ending pos.
        z_marker_weight_sets[s] = to_row_vector(marker_weight_prior_effect[start_pos:end_pos]);
      }
    }

    if (assoc_cv_total == 1 && marker_weight_set_cv_total > 0) {
      marker_weights_eff_cv_total = marker_weights_cv_total + to_vector(z_marker_weight_sets[marker_weight_set_cv_total]);
    }
    if (assoc_cs_total == 1 && marker_weight_set_cs_total > 0) {
      marker_weights_eff_cs_total = marker_weights_cs_total + to_vector(z_marker_weight_sets[marker_weight_set_cs_total]);
    }
    if (assoc_cv_marker == 1 && marker_weight_set_cv_marker > 0) {
      marker_weights_eff_cv_marker = marker_weights_cv_marker + to_vector(z_marker_weight_sets[marker_weight_set_cv_marker]);
    }
    if (assoc_cs_marker == 1 && marker_weight_set_cs_marker > 0) {
      marker_weights_eff_cs_marker = marker_weights_cs_marker + to_vector(z_marker_weight_sets[marker_weight_set_cs_marker]);
    }
  }

  /* -------------------- effective spline coefficients for transforms */
  // For every Stan-estimated monotone transform, we estimate only its SHAPE.
  // The first ordinate is anchored at zero and the absolute span is one. An
  // increasing I-spline or piecewise-linear curve ends at one; a decreasing
  // piecewise-linear curve ends at minus one. Decreasing I-spline modes retain
  // increasing coefficients because their evaluator reflects the basis. This
  // avoids confounding transform scale with the association coefficient.
  // Implementation detail:
  // - piecewise-linear modes sample positive increments directly as a simplex,
  // - established I-spline modes retain softmax-transformed latent increments,
  // - in either case the increments sum to 1 exactly,
  // - cumulative sums therefore give a stable unit-span monotone curve
  //   without any divide-by-a-nearly-zero normalisation step.
  vector[n_coeff_cv] coeff_cv_eff = coeff_cv; // Role: coefficients current value eff.
  vector[n_coeff_cs] coeff_cs_eff = coeff_cs; // Role: coefficients current slope eff.
  matrix[M_corr, n_coeff_corr] coeff_corr_eff = coeff_corr; // Role: coefficients correlation eff.
  matrix[M_vcov, n_coeff_vcov] coeff_vcov_eff = coeff_vcov; // Role: coefficients covariance eff.
  vector[n_coeff_cv_mean] coeff_cv_mean_eff = coeff_cv_mean; // Role: coefficients current value mean eff.
  vector[n_coeff_cv_marker] coeff_cv_marker_eff = coeff_cv_marker; // Role: coefficients current value marker eff.
  vector[n_coeff_cs_mean] coeff_cs_mean_eff = coeff_cs_mean; // Role: coefficients current slope mean eff.
  vector[n_coeff_cs_marker] coeff_cs_marker_eff = coeff_cs_marker; // Role: coefficients current slope marker eff.

  if (estimate_spline_cv == 1 && n_coeff_cv > 1 && n_free_spline_cv == n_coeff_cv - 1) {
    vector[n_coeff_cv - 1] delta; // Role: increment.
    if (tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7)
      delta = pwlin_simplex_cv;
    else
      delta = softmax(z_spline_cv);
    coeff_cv_eff[1] = 0;
    for (j in 2:n_coeff_cv)
      coeff_cv_eff[j] = coeff_cv_eff[j - 1] + ((tf_mode_cv_tot == 7) ? -delta[j - 1] : delta[j - 1]);
  }
  if (estimate_spline_cs == 1 && n_coeff_cs > 1 && n_free_spline_cs == n_coeff_cs - 1) {
    vector[n_coeff_cs - 1] delta; // Role: increment.
    if (tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7)
      delta = pwlin_simplex_cs;
    else
      delta = softmax(z_spline_cs);
    coeff_cs_eff[1] = 0;
    for (j in 2:n_coeff_cs)
      coeff_cs_eff[j] = coeff_cs_eff[j - 1] + ((tf_mode_cs_tot == 7) ? -delta[j - 1] : delta[j - 1]);
  }
  if (estimate_spline_corr == 1 && n_coeff_corr > 1 && n_free_spline_corr == n_coeff_corr - 1) {
    for (m in 1:M_corr) {
      vector[n_coeff_corr - 1] delta; // Role: increment.
      if (tf_mode_corr == 3 || tf_mode_corr == 7)
        delta = pwlin_simplex_corr[m];
      else
        delta = softmax(to_vector(row(z_spline_corr, m)));
      coeff_corr_eff[m, 1] = 0;
      for (j in 2:n_coeff_corr)
        coeff_corr_eff[m, j] = coeff_corr_eff[m, j - 1] + ((tf_mode_corr == 7) ? -delta[j - 1] : delta[j - 1]);
    }
  }
  if (estimate_spline_vcov == 1 && n_coeff_vcov > 1 && n_free_spline_vcov == n_coeff_vcov - 1) {
    for (m in 1:M_vcov) {
      vector[n_coeff_vcov - 1] delta; // Role: increment.
      if (tf_mode_vcov == 3 || tf_mode_vcov == 7)
        delta = pwlin_simplex_vcov[m];
      else
        delta = softmax(to_vector(row(z_spline_vcov, m)));
      coeff_vcov_eff[m, 1] = 0;
      for (j in 2:n_coeff_vcov)
        coeff_vcov_eff[m, j] = coeff_vcov_eff[m, j - 1] + ((tf_mode_vcov == 7) ? -delta[j - 1] : delta[j - 1]);
    }
  }
  if (estimate_spline_cv_mean == 1 && n_coeff_cv_mean > 1 && n_free_spline_cv_mean == n_coeff_cv_mean - 1) {
    vector[n_coeff_cv_mean - 1] delta; // Role: increment.
    if (tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7)
      delta = pwlin_simplex_cv_mean;
    else
      delta = softmax(z_spline_cv_mean);
    coeff_cv_mean_eff[1] = 0;
    for (j in 2:n_coeff_cv_mean)
      coeff_cv_mean_eff[j] = coeff_cv_mean_eff[j - 1] + ((tf_mode_cv_mean == 7) ? -delta[j - 1] : delta[j - 1]);
  }
  if (estimate_spline_cv_marker == 1 && n_coeff_cv_marker > 1 && n_free_spline_cv_marker == n_coeff_cv_marker - 1) {
    vector[n_coeff_cv_marker - 1] delta; // Role: increment.
    if (tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7)
      delta = pwlin_simplex_cv_marker;
    else
      delta = softmax(z_spline_cv_marker);
    coeff_cv_marker_eff[1] = 0;
    for (j in 2:n_coeff_cv_marker)
      coeff_cv_marker_eff[j] = coeff_cv_marker_eff[j - 1] + ((tf_mode_cv_marker == 7) ? -delta[j - 1] : delta[j - 1]);
  }
  if (estimate_spline_cs_mean == 1 && n_coeff_cs_mean > 1 && n_free_spline_cs_mean == n_coeff_cs_mean - 1) {
    vector[n_coeff_cs_mean - 1] delta; // Role: increment.
    if (tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7)
      delta = pwlin_simplex_cs_mean;
    else
      delta = softmax(z_spline_cs_mean);
    coeff_cs_mean_eff[1] = 0;
    for (j in 2:n_coeff_cs_mean)
      coeff_cs_mean_eff[j] = coeff_cs_mean_eff[j - 1] + ((tf_mode_cs_mean == 7) ? -delta[j - 1] : delta[j - 1]);
  }
  if (estimate_spline_cs_marker == 1 && n_coeff_cs_marker > 1 && n_free_spline_cs_marker == n_coeff_cs_marker - 1) {
    vector[n_coeff_cs_marker - 1] delta; // Role: increment.
    if (tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7)
      delta = pwlin_simplex_cs_marker;
    else
      delta = softmax(z_spline_cs_marker);
    coeff_cs_marker_eff[1] = 0;
    for (j in 2:n_coeff_cs_marker)
      coeff_cs_marker_eff[j] = coeff_cs_marker_eff[j - 1] + ((tf_mode_cs_marker == 7) ? -delta[j - 1] : delta[j - 1]);
  }

  /* -------------------- fit-only affine shift for functional transforms */
  vector[estimate_iota_intercept_cv] iota_intercept_cv_eff = rep_vector(0, estimate_iota_intercept_cv); // Role: affine transformation intercept current value eff.
  vector[estimate_iota_slope_cv] iota_slope_cv_eff = rep_vector(1, estimate_iota_slope_cv); // Role: affine transformation slope current value eff.
  vector[estimate_iota_intercept_cs] iota_intercept_cs_eff = rep_vector(0, estimate_iota_intercept_cs); // Role: affine transformation intercept current slope eff.
  vector[estimate_iota_slope_cs] iota_slope_cs_eff = rep_vector(1, estimate_iota_slope_cs); // Role: affine transformation slope current slope eff.
  vector[M_corr * estimate_iota_intercept_corr] iota_intercept_corr_eff = rep_vector(0, M_corr * estimate_iota_intercept_corr); // Role: affine transformation intercept correlation eff.
  vector[M_corr * estimate_iota_slope_corr] iota_slope_corr_eff = rep_vector(1, M_corr * estimate_iota_slope_corr); // Role: affine transformation slope correlation eff.
  vector[M_vcov * estimate_iota_intercept_vcov] iota_intercept_vcov_eff = rep_vector(0, M_vcov * estimate_iota_intercept_vcov); // Role: affine transformation intercept covariance eff.
  vector[M_vcov * estimate_iota_slope_vcov] iota_slope_vcov_eff = rep_vector(1, M_vcov * estimate_iota_slope_vcov); // Role: affine transformation slope covariance eff.
  vector[estimate_iota_intercept_cv_mean] iota_intercept_cv_mean_eff = rep_vector(0, estimate_iota_intercept_cv_mean); // Role: affine transformation intercept current value mean eff.
  vector[estimate_iota_slope_cv_mean] iota_slope_cv_mean_eff = rep_vector(1, estimate_iota_slope_cv_mean); // Role: affine transformation slope current value mean eff.
  vector[estimate_iota_intercept_cv_marker] iota_intercept_cv_marker_eff = rep_vector(0, estimate_iota_intercept_cv_marker); // Role: affine transformation intercept current value marker eff.
  vector[estimate_iota_slope_cv_marker] iota_slope_cv_marker_eff = rep_vector(1, estimate_iota_slope_cv_marker); // Role: affine transformation slope current value marker eff.
  vector[estimate_iota_intercept_cs_mean] iota_intercept_cs_mean_eff = rep_vector(0, estimate_iota_intercept_cs_mean); // Role: affine transformation intercept current slope mean eff.
  vector[estimate_iota_slope_cs_mean] iota_slope_cs_mean_eff = rep_vector(1, estimate_iota_slope_cs_mean); // Role: affine transformation slope current slope mean eff.
  vector[estimate_iota_intercept_cs_marker] iota_intercept_cs_marker_eff = rep_vector(0, estimate_iota_intercept_cs_marker); // Role: affine transformation intercept current slope marker eff.
  vector[estimate_iota_slope_cs_marker] iota_slope_cs_marker_eff = rep_vector(1, estimate_iota_slope_cs_marker); // Role: affine transformation slope current slope marker eff.

  {
    vector[n_iota_prior] iota_raw; // all affine-shift raw parameters in their documented block order
    vector[n_iota_prior] iota_prior_effect; // scientific affine shifts after the selected prior transformation
    int iota_position = 1; // next free position in the packed affine-shift vector

    if (estimate_iota_intercept_cv > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cv - 1] = z_iota_intercept_cv; iota_position += estimate_iota_intercept_cv; }
    if (estimate_iota_slope_cv > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cv - 1] = z_iota_slope_cv; iota_position += estimate_iota_slope_cv; }
    if (estimate_iota_intercept_cs > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cs - 1] = z_iota_intercept_cs; iota_position += estimate_iota_intercept_cs; }
    if (estimate_iota_slope_cs > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cs - 1] = z_iota_slope_cs; iota_position += estimate_iota_slope_cs; }
    if (M_corr * estimate_iota_intercept_corr > 0) { iota_raw[iota_position : iota_position + M_corr * estimate_iota_intercept_corr - 1] = z_iota_intercept_corr; iota_position += M_corr * estimate_iota_intercept_corr; }
    if (M_corr * estimate_iota_slope_corr > 0) { iota_raw[iota_position : iota_position + M_corr * estimate_iota_slope_corr - 1] = z_iota_slope_corr; iota_position += M_corr * estimate_iota_slope_corr; }
    if (M_vcov * estimate_iota_intercept_vcov > 0) { iota_raw[iota_position : iota_position + M_vcov * estimate_iota_intercept_vcov - 1] = z_iota_intercept_vcov; iota_position += M_vcov * estimate_iota_intercept_vcov; }
    if (M_vcov * estimate_iota_slope_vcov > 0) { iota_raw[iota_position : iota_position + M_vcov * estimate_iota_slope_vcov - 1] = z_iota_slope_vcov; iota_position += M_vcov * estimate_iota_slope_vcov; }
    if (estimate_iota_intercept_cv_mean > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cv_mean - 1] = z_iota_intercept_cv_mean; iota_position += estimate_iota_intercept_cv_mean; }
    if (estimate_iota_slope_cv_mean > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cv_mean - 1] = z_iota_slope_cv_mean; iota_position += estimate_iota_slope_cv_mean; }
    if (estimate_iota_intercept_cv_marker > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cv_marker - 1] = z_iota_intercept_cv_marker; iota_position += estimate_iota_intercept_cv_marker; }
    if (estimate_iota_slope_cv_marker > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cv_marker - 1] = z_iota_slope_cv_marker; iota_position += estimate_iota_slope_cv_marker; }
    if (estimate_iota_intercept_cs_mean > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cs_mean - 1] = z_iota_intercept_cs_mean; iota_position += estimate_iota_intercept_cs_mean; }
    if (estimate_iota_slope_cs_mean > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cs_mean - 1] = z_iota_slope_cs_mean; iota_position += estimate_iota_slope_cs_mean; }
    if (estimate_iota_intercept_cs_marker > 0) { iota_raw[iota_position : iota_position + estimate_iota_intercept_cs_marker - 1] = z_iota_intercept_cs_marker; iota_position += estimate_iota_intercept_cs_marker; }
    if (estimate_iota_slope_cs_marker > 0) { iota_raw[iota_position : iota_position + estimate_iota_slope_cs_marker - 1] = z_iota_slope_cs_marker; }

    iota_prior_effect = joinme_prior_transform(
      iota_raw,
      prior_iota_family,
      prior_iota_mu,
      prior_iota_scale,
      horseshoe_local_iota,
      horseshoe_global_iota,
      horseshoe_slab_iota,
      prior_iota_slab_scale
    );

    iota_position = 1;
    if (estimate_iota_intercept_cv > 0) { iota_intercept_cv_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cv - 1]; iota_position += estimate_iota_intercept_cv; }
    if (estimate_iota_slope_cv > 0) { iota_slope_cv_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cv - 1]; iota_position += estimate_iota_slope_cv; }
    if (estimate_iota_intercept_cs > 0) { iota_intercept_cs_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cs - 1]; iota_position += estimate_iota_intercept_cs; }
    if (estimate_iota_slope_cs > 0) { iota_slope_cs_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cs - 1]; iota_position += estimate_iota_slope_cs; }
    if (M_corr * estimate_iota_intercept_corr > 0) { iota_intercept_corr_eff = iota_prior_effect[iota_position : iota_position + M_corr * estimate_iota_intercept_corr - 1]; iota_position += M_corr * estimate_iota_intercept_corr; }
    if (M_corr * estimate_iota_slope_corr > 0) { iota_slope_corr_eff = iota_prior_effect[iota_position : iota_position + M_corr * estimate_iota_slope_corr - 1]; iota_position += M_corr * estimate_iota_slope_corr; }
    if (M_vcov * estimate_iota_intercept_vcov > 0) { iota_intercept_vcov_eff = iota_prior_effect[iota_position : iota_position + M_vcov * estimate_iota_intercept_vcov - 1]; iota_position += M_vcov * estimate_iota_intercept_vcov; }
    if (M_vcov * estimate_iota_slope_vcov > 0) { iota_slope_vcov_eff = iota_prior_effect[iota_position : iota_position + M_vcov * estimate_iota_slope_vcov - 1]; iota_position += M_vcov * estimate_iota_slope_vcov; }
    if (estimate_iota_intercept_cv_mean > 0) { iota_intercept_cv_mean_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cv_mean - 1]; iota_position += estimate_iota_intercept_cv_mean; }
    if (estimate_iota_slope_cv_mean > 0) { iota_slope_cv_mean_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cv_mean - 1]; iota_position += estimate_iota_slope_cv_mean; }
    if (estimate_iota_intercept_cv_marker > 0) { iota_intercept_cv_marker_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cv_marker - 1]; iota_position += estimate_iota_intercept_cv_marker; }
    if (estimate_iota_slope_cv_marker > 0) { iota_slope_cv_marker_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cv_marker - 1]; iota_position += estimate_iota_slope_cv_marker; }
    if (estimate_iota_intercept_cs_mean > 0) { iota_intercept_cs_mean_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cs_mean - 1]; iota_position += estimate_iota_intercept_cs_mean; }
    if (estimate_iota_slope_cs_mean > 0) { iota_slope_cs_mean_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cs_mean - 1]; iota_position += estimate_iota_slope_cs_mean; }
    if (estimate_iota_intercept_cs_marker > 0) { iota_intercept_cs_marker_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_intercept_cs_marker - 1]; iota_position += estimate_iota_intercept_cs_marker; }
    if (estimate_iota_slope_cs_marker > 0) { iota_slope_cs_marker_eff = iota_prior_effect[iota_position : iota_position + estimate_iota_slope_cs_marker - 1]; }
  }

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
    wbar_i[i] = diag_pre_multiply(row_scale_idm, L_i[i]) * zbar[i];
  
  /* -------------------- Effective association coefficients (flags applied) */
  vector[6 + M_corr + M_vcov] alpha_prior_raw; // six scalar association latents followed by corr and vcov vectors
  vector[6 + M_corr + M_vcov] alpha_prior_effect; // association latents after their block-specific prior transformation
  alpha_prior_raw[1] = z_alpha_cv_total;
  alpha_prior_raw[2] = z_alpha_cs_total;
  alpha_prior_raw[3] = z_alpha_cv_mean;
  alpha_prior_raw[4] = z_alpha_cs_mean;
  alpha_prior_raw[5] = z_alpha_cv_marker;
  alpha_prior_raw[6] = z_alpha_cs_marker;
  if (M_corr > 0) alpha_prior_raw[7 : 6 + M_corr] = z_alpha_corr;
  if (M_vcov > 0) alpha_prior_raw[7 + M_corr : 6 + M_corr + M_vcov] = z_alpha_vcov;
  alpha_prior_effect = joinme_prior_transform(
    alpha_prior_raw,
    prior_alpha_family,
    prior_alpha_mu,
    prior_alpha_scale,
    horseshoe_local_alpha,
    horseshoe_global_alpha,
    horseshoe_slab_alpha,
    prior_alpha_slab_scale
  );

  // Non-centred association construction with flexible sign constraints:
  // - when weighted marker-based association terms share one marker-weight
  //   structure, only the first active weighted term is constrained positive;
  // - when weighted marker-based association terms use separate weight
  //   structures, every active weighted term is constrained positive.
  real z_alpha_cv_total_constrained = alpha_prior_effect[1]; // transformed latent for total current-value association
  real z_alpha_cs_total_constrained = alpha_prior_effect[2]; // transformed latent for total current-slope association
  real z_alpha_cv_marker_constrained = alpha_prior_effect[5]; // transformed latent for marker current-value association
  real z_alpha_cs_marker_constrained = alpha_prior_effect[6]; // transformed latent for marker current-slope association

  if ((assoc_cv_total + assoc_cs_total + assoc_cv_marker + assoc_cs_marker) > 0) {
    if (shared_marker_weights == 1) {
      if (assoc_cv_total == 1) {
        z_alpha_cv_total_constrained = abs(alpha_prior_effect[1]);
      } else if (assoc_cs_total == 1) {
        z_alpha_cs_total_constrained = abs(alpha_prior_effect[2]);
      } else if (assoc_cv_marker == 1) {
        z_alpha_cv_marker_constrained = abs(alpha_prior_effect[5]);
      } else if (assoc_cs_marker == 1) {
        z_alpha_cs_marker_constrained = abs(alpha_prior_effect[6]);
      }
    } else {
      if (assoc_cv_total == 1) z_alpha_cv_total_constrained = abs(alpha_prior_effect[1]);
      if (assoc_cs_total == 1) z_alpha_cs_total_constrained = abs(alpha_prior_effect[2]);
      if (assoc_cv_marker == 1) z_alpha_cv_marker_constrained = abs(alpha_prior_effect[5]);
      if (assoc_cs_marker == 1) z_alpha_cs_marker_constrained = abs(alpha_prior_effect[6]);
    }
  }

  /*
   * These quantities are the association coefficients used by the hazard.
   * They are not multiplied by an additional unidentified random scale: the
   * location and scale supplied in `jm_prior(alpha = ...)` therefore describe
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
