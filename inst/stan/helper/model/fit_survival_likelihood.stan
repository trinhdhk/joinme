  /**
   * @brief Survival likelihood with association terms and Gaussian-Kronrod integration.
   *
   * For each subject i:
   * 1. Compute baseline log-hazard at 15 GK nodes and at event time
   * 2. Build raw association components (current value, slope, variance covariance)
   * 3. Apply transformations to association features according to mode
   * 4. Accumulate event term if subject had event (d_event[i] == 1)
   * 5. Subtract cumulative hazard to complete log-likelihood
   *
  * Association features are computed using:
  * - Current value total (mean + marker effects; marker averages are weighted)
  * - Current slope total (mean slope + marker slopes; same weights apply)
  * - Variance covariance feature (id-specific variance structure)
   * Each can have its own transformation (functional, spline, etc.)
   */

  for (i in 1 : n_id) {
    // We build association terms once per subject, then reuse them for each cause.
    
    // raw components at GK nodes
    vector[15] cvm_now;
    vector[15] cvm_fwd;
    vector[15] cvk_now;
    vector[15] cvk_fwd;
    
    for (j in 1 : 15) {
      // mean component: X + id RE
      cvm_now[j] = dot_product(X_gk_now[i][j], beta_scaled)
                   + dot_product(Z_id_gk_now[i][j], u_id[i]);
      cvm_fwd[j] = dot_product(X_gk_fwd[i][j], beta_scaled)
                   + dot_product(Z_id_gk_fwd[i][j], u_id[i]);
      
      // marker-average component: optional marker-only + marker-by-id
      // Note: vbar and wbar_i are weighted means across markers using marker_weights.
      // These weighted averages feed both current value (CV) and current slope (CS).
      real mk_part_now = 0;
      real mk_part_fwd = 0;
      if (R_mk > 0) {
        mk_part_now = dot_product(Z_mk_gk_now[i][j], vbar);
        mk_part_fwd = dot_product(Z_mk_gk_fwd[i][j], vbar);
      }
      
      cvk_now[j] = mk_part_now + dot_product(Z_idm_gk_now[i][j], wbar_i[i]);
      cvk_fwd[j] = mk_part_fwd + dot_product(Z_idm_gk_fwd[i][j], wbar_i[i]);
    }
    
    // FD slopes (marker weights are already folded into cvk via vbar/wbar_i):
    // eta_fd returns slope with respect to the scaled-time domain used in the GK designs.
    // Convert to original-time slope by dividing by tmax.
    vector[15] csm_raw = eta_fd(cvm_now, cvm_fwd, eps_fd) / tmax;
    vector[15] csk_raw = eta_fd(cvk_now, cvk_fwd, eps_fd) / tmax;
    
    // vcov feature
    vector[15] vcov_raw = eta_vcov_varonly_weighted_const(15, L_i[i],
                            a_vcov_var);
    
    // -------------------- Association transforms
    vector[15] cv_tot = cvm_now + cvk_now;
    vector[15] cs_tot = csm_raw + csk_raw;
    
    vector[15] cv_tot_tf = apply_transform_vector(
      cv_tot,
      tf_mode_cv_tot,
      functional_ops_cv,
      const_data_cv,
      knots_cv,
      coeff_cv,
      spline_degree_cv
    );
    vector[15] cv_mean_tf = apply_transform_vector(
      cvm_now,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
      coeff_cv_mean,
      spline_degree_cv_mean
    );
    vector[15] cv_marker_tf = apply_transform_vector(
      cvk_now,
      tf_mode_cv_marker,
      functional_ops_cv_marker,
      const_data_cv_marker,
      knots_cv_marker,
      coeff_cv_marker,
      spline_degree_cv_marker
    );
    
    vector[15] cs_tot_tf = apply_transform_vector(
      cs_tot,
      tf_mode_cs_tot,
      functional_ops_cs,
      const_data_cs,
      knots_cs,
      coeff_cs,
      spline_degree_cs
    );
    vector[15] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
      coeff_cs_mean,
      spline_degree_cs_mean
    );
    vector[15] cs_marker_tf = apply_transform_vector(
      csk_raw,
      tf_mode_cs_marker,
      functional_ops_cs_marker,
      const_data_cs_marker,
      knots_cs_marker,
      coeff_cs_marker,
      spline_degree_cs_marker
    );
    
    vector[15] vcov_tf = apply_transform_vector(
      vcov_raw,
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov,
      spline_degree_vcov
    );
    
    vector[15] eta_assoc_nodes = a_cv_total * cv_tot_tf
                           + a_cv_mean * cv_mean_tf
                           + a_cv_marker * cv_marker_tf
                           + a_cs_total * cs_tot_tf
                           + a_cs_mean * cs_mean_tf
                           + a_cs_marker * cs_marker_tf
                           + vcov_tf;
    
    // event term (use the observed event type)
    if (d_event[i] == 1) {
      int k_ev = event_type[i];
      real eta_w_ev = 0;
      if (p_w > 0)
        eta_w_ev = dot_product(to_vector(W[i]'), gamma_w[k_ev]);

      real log_h0_S = log_h0_intercept[k_ev]
                      + dot_product(Bs_event_c[i], bs_gamma_c[k_ev]);
      
      real cvm_S = dot_product(X_event_now[i], beta_scaled)
                   + dot_product(Z_id_event_now[i], u_id[i]);
      real cvm_S_fwd = dot_product(X_event_fwd[i], beta_scaled)
                       + dot_product(Z_id_event_fwd[i], u_id[i]);
      
      real mk_part_S = 0;
      real mk_part_S_fwd = 0;
      if (R_mk > 0) {
        mk_part_S = dot_product(Z_mk_event_now[i], vbar);
        mk_part_S_fwd = dot_product(Z_mk_event_fwd[i], vbar);
      }
      
      real cvk_S = mk_part_S + dot_product(Z_idm_event_now[i], wbar_i[i]);
      real cvk_S_fwd = mk_part_S_fwd
                       + dot_product(Z_idm_event_fwd[i], wbar_i[i]);
      
      // slopes at event time (convert to original-time slope)
      real csm_S_raw = ((cvm_S_fwd - cvm_S) / eps_fd) / tmax;
      real csk_S_raw = ((cvk_S_fwd - cvk_S) / eps_fd) / tmax;
      
      // vcov scalar at event time
      real vcov_S_raw = 0;
      for (q in 1 : Q_idm) {
        real var_q = 0;
        for (k in 1 : q) 
          var_q += square(L_i[i][q, k]);
        vcov_S_raw += a_vcov_var[q] * log(var_q);
      }
      
      // totals + transforms
      real cv_S_tot = cvm_S + cvk_S;
      real cs_S_tot = csm_S_raw + csk_S_raw;
      
      real cv_S_tot_tf = apply_transform_scalar(
        cv_S_tot,
        tf_mode_cv_tot,
        functional_ops_cv,
        const_data_cv,
        knots_cv,
        coeff_cv,
        spline_degree_cv
      );
      real cv_S_mean_tf = apply_transform_scalar(
        cvm_S,
        tf_mode_cv_mean,
        functional_ops_cv_mean,
        const_data_cv_mean,
        knots_cv_mean,
        coeff_cv_mean,
        spline_degree_cv_mean
      );
      real cv_S_marker_tf = apply_transform_scalar(
        cvk_S,
        tf_mode_cv_marker,
        functional_ops_cv_marker,
        const_data_cv_marker,
        knots_cv_marker,
        coeff_cv_marker,
        spline_degree_cv_marker
      );
      
      real cs_S_tot_tf = apply_transform_scalar(
        cs_S_tot,
        tf_mode_cs_tot,
        functional_ops_cs,
        const_data_cs,
        knots_cs,
        coeff_cs,
        spline_degree_cs
      );
      real cs_S_mean_tf = apply_transform_scalar(
        csm_S_raw,
        tf_mode_cs_mean,
        functional_ops_cs_mean,
        const_data_cs_mean,
        knots_cs_mean,
        coeff_cs_mean,
        spline_degree_cs_mean
      );
      real cs_S_marker_tf = apply_transform_scalar(
        csk_S_raw,
        tf_mode_cs_marker,
        functional_ops_cs_marker,
        const_data_cs_marker,
        knots_cs_marker,
        coeff_cs_marker,
        spline_degree_cs_marker
      );
      
      real vcov_S_tf = apply_transform_scalar(
        vcov_S_raw,
        tf_mode_vcov,
        functional_ops_vcov,
        const_data_vcov,
        knots_vcov,
        coeff_vcov,
        spline_degree_vcov
      );
      
      real eta_assoc_S = a_cv_total * cv_S_tot_tf
                   + a_cv_mean * cv_S_mean_tf
                   + a_cv_marker * cv_S_marker_tf
                   + a_cs_total * cs_S_tot_tf
                   + a_cs_mean * cs_S_mean_tf
                   + a_cs_marker * cs_S_marker_tf
                   + vcov_S_tf;
      
      target += log_h0_S + eta_w_ev + eta_assoc_S;
    }

    // cumulative hazard: sum of cause-specific hazards
    vector[15] log_h_total;
    for (j in 1 : 15) {
      vector[K_event] log_h_cause;
      for (k_ev in 1 : K_event) {
        real eta_w_k = 0;
        if (p_w > 0)
          eta_w_k = dot_product(to_vector(W[i]'), gamma_w[k_ev]);
        log_h_cause[k_ev] = log_h0_intercept[k_ev]
                            + dot_product(Bs_gk_c[i][j], bs_gamma_c[k_ev])
                            + eta_w_k
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    target += -cumhaz(S_event[i], log_h_total, rep_vector(0.0, 15));
  }
