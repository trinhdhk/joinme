/**
 * @file longitudinal/transformed_parameters/fit.stan
 * @brief Construct longitudinal effects and subject-specific covariance factors.
 *
 * @details
 * Population and random-effect parameters already use the original study-time
 * basis supplied by R. The fragment constructs subject effects, marker effects,
 * subject-by-marker latent effects, subject-specific Cholesky factors and their
 * realised random effects without any coefficient rescaling.
 */
  /* -------------------- id RE covariance */
  matrix[R_id, R_id] L_u; // Cholesky factor for subject effects on original time
  if (indep_id_re == 1) 
    L_u = diag_matrix(tau_u);
  else 
    L_u = diag_pre_multiply(tau_u, Lcorr_u);
  
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
      rep_vector(0, D * R_mk),
      rep_vector(1, D * R_mk),
      horseshoe_local_marker,
      horseshoe_global_marker,
      horseshoe_slab_marker,
      prior_marker_slab_scale
    );
    if (indep_marker_re == 1) 
      L_v = diag_matrix(tau_v);
    else 
      L_v = diag_pre_multiply(tau_v, Lcorr_v);
    
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
  
  /* -------------------- subject-specific Cholesky factors L_i = SD_i * K_i */
  array[n_id] matrix[Q_idm, Q_idm] L_i; // subject-specific covariance Cholesky factor
  for (i in 1 : n_id) {
    matrix[Q_idm, Q_idm] Ki = rep_matrix(0.0, Q_idm, Q_idm); // row-normalised Cholesky correlation factor
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

    // Step 2: rebuild K_i under the independent SD+correlation parameterisation.
    //
    // Each tanh-linked predictor is a partial correlation. Multiplying by the
    // square root of the variance left by preceding coordinates gives a row
    // of a valid Cholesky correlation factor. Every row of K_i therefore has
    // norm one, so sd_i[r] remains the marginal standard deviation of random
    // effect r rather than changing with its correlation predictors. The
    // covariance factor is L_i = diag_pre_multiply(sd_i, K_i).
    //
    // For an off-diagonal coordinate, lambda_L[m] is the residual scale before
    // tanh; it is therefore not itself a correlation standard deviation.
    for (r in 1 : Q_idm) {
      if (r == 1) {
        Ki[1, 1] = 1.0;
        if (indep_idmarker_cov == 0)
          m_pos += 1;
      } else if (indep_idmarker_cov == 1) {
        Ki[r, r] = 1.0;
        m_pos += 1;
      } else {
        real remaining_scale = 1.0; // square-root scale not assigned to preceding coordinates in this row
        for (c in 1 : r) {
          if (c < r) {
            int correlation_coordinate = ((r - 1) * (r - 2)) %/% 2 + c; // row-major off-diagonal index: (2,1), (3,1), (3,2), ...
            real lp = alpha_L[m_pos] + dot_product(beta_L_corr[correlation_coordinate], to_vector(Xcov_corr[i]'))
              + lambda_L[m_pos] * z_L[i][m_pos];
            real partial_correlation = tanh(lp); // partial correlation for coordinate (r,c)
            Ki[r, c] = remaining_scale * partial_correlation;
            remaining_scale *= sqrt(fmax(1e-12, 1.0 - square(partial_correlation)));
          } else {
            Ki[r, r] = remaining_scale;
          }
          m_pos += 1;
        }
      }
    }

    Li = diag_pre_multiply(sd_i, Ki);

    L_i[i] = Li;
  }
  
  /* -------------------- marker-by-subject effects: w_idm[i,d] = L_i[i] * z_w[i,d] */
  array[n_id, D] vector[Q_idm] w_idm; // realised marker-by-subject effects on original time
  for (i in 1 : n_id) {
    for (d in 1 : D)
      w_idm[i, d] = L_i[i] * z_w[i, d];
  }
  
