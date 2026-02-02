/**
 * @file joinme_dynpred.stan
 * @brief Dynamic prediction for joint longitudinal + survival models.
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/functions/functional_transform.stanfunctions
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions
}

data {
  #include helper/data/joinme_dynpred_common.stan

  int flag_indep_id_re;
  int flag_indep_marker_re;
  int flag_indep_marker_byid_latent_re;
  int flag_indep_idmarker_cov;
  int flag_allow_marker_crosscorr;

  // Mappings for covariance regression (pre-calculated indices)
  array[num_unique_cov_entries] int idx_row_cov;
  array[num_unique_cov_entries] int idx_col_cov;
}

transformed data {
}
parameters {
  #include helper/parameters/joinme_dynpred_common.stan
}
model {
  // Posterior Sampling of Random Effects for the New Subject
  // We re-use the likelihood logic from generated quantities (duplicated for validity in 'model' block)
  // Logic is minimized here to just what is needed for sampling REs.
  
  for (i in 1 : n_draws) {
    // Priors
    z_u[i] ~ std_normal();
    if (n_random_marker > 0) {
      for (d in 1 : n_marker_types) 
        z_v[i, d] ~ std_normal();
    }
    for (d in 1 : n_marker_types) 
      z_w_lat[i, d] ~ std_normal();
    z_L[i] ~ std_normal();
    
    // ------------------------------------------
    // Construct REs
    // ------------------------------------------
    matrix[n_random_id, n_random_id] L_u;
    if (flag_indep_id_re == 1) 
      L_u = diag_matrix(tau_id[i]);
    else 
      L_u = diag_pre_multiply(tau_id[i], Lcorr_id[i]);
    vector[n_random_id] u_id = L_u * z_u[i];
    
    matrix[n_random_marker, n_random_marker] L_v;
    if (n_random_marker > 0) {
      if (flag_indep_marker_re == 1) 
        L_v = diag_matrix(tau_marker[i]);
      else 
        L_v = diag_pre_multiply(tau_marker[i], Lcorr_marker[i]);
    }
    array[n_marker_types] vector[n_random_marker] v_marker;
    if (n_random_marker > 0) {
      for (d in 1 : n_marker_types) 
        v_marker[d] = L_v * z_v[i, d];
    }
    
    matrix[n_random_marker_id, n_random_marker_id] L_w;
    if (flag_indep_marker_byid_latent_re == 1) 
      L_w = diag_matrix(tau_marker_id[i]);
    else 
      L_w = diag_pre_multiply(tau_marker_id[i], Lcorr_marker_id[i]);
    array[n_marker_types] vector[n_random_marker_id] z_w;
    for (d in 1 : n_marker_types) {
      vector[n_random_marker_id] cross = rep_vector(0.0, n_random_marker_id);
      if (flag_allow_marker_crosscorr == 1 && n_random_marker > 0) 
        cross = B_cross[i] * v_marker[d];
      z_w[d] = cross + L_w * z_w_lat[i, d];
    }
    
    real u_Lk = tau_vcov_reg[i] * z_L[i];
    matrix[n_random_marker_id, n_random_marker_id] Li = rep_matrix(0.0,
                                                                   n_random_marker_id,
                                                                   n_random_marker_id);
    for (m in 1 : num_unique_cov_entries) {
      int r_ = idx_row_cov[m];
      int c_ = idx_col_cov[m];
      vector[n_cov_vcov] bL_m = beta_vcov_reg_flat[i][((m - 1) * n_cov_vcov
                                                       + 1) : (m * n_cov_vcov)];
      real lp = alpha_vcov_reg[i][m] + dot_product(bL_m, vec_cov_vcov)
                + lambda_vcov_reg[i][m] * u_Lk;
      Li[r_, c_] = (r_ == c_) ? log1p_exp(lp) : lp;
    }
    
    array[n_marker_types] vector[n_random_marker_id] w_idscaled;
    for (d in 1 : n_marker_types) 
      w_idscaled[d] = Li * z_w[d];
    
    // ------------------------------------------
    // Longitudinal Likelihood
    // ------------------------------------------
    for (n in 1 : n_obs_long) {
      int d = idx_marker_obs[n];
      real eta = dot_product(mat_fixed_obs[n], beta_fixed[i])
                 + dot_product(mat_id_obs[n], u_id)
                 + ((n_random_marker > 0)
                    ? dot_product(mat_marker_obs[n], v_marker[d]) : 0.0)
                 + dot_product(mat_marker_id_obs[n], w_idscaled[d]);
      
      if (family_long[d] == 1) {
        real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[i]))
                   : ((flag_resid_dim == 1) ? sigma_marker_specific[i][d]
                      : sigma_y_shared[i]);
        target += normal_lpdf(y_real[n] | eta, sig);
      } else if (family_long[d] == 2) {
        real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[i]))
                   : ((flag_resid_dim == 1) ? sigma_marker_specific[i][d]
                      : sigma_y_shared[i]);
        real nu = (P_nu > 0) ? (2 + exp(dot_product(X_nu_obs[n], beta_nu[i])))
                  : nu_marker[i][d];
        target += student_t_lpdf(y_real[n] | nu, eta, sig);
      } else if (family_long[d] == 3) {
        target += bernoulli_logit_lpmf(y_int[n] | eta);
      } else if (family_long[d] == 4) {
        target += binomial_logit_lpmf(y_int[n] | trials_obs[n], eta);
      } else if (family_long[d] == 5) {
        target += poisson_log_lpmf(y_int[n] | eta);
      } else if (family_long[d] == 6) {
        real phi = (P_phi > 0) ? exp(dot_product(X_phi_obs[n], beta_phi[i]))
                   : phi_nb_marker[i][d];
        target += neg_binomial_2_log_lpmf(y_int[n] | eta, phi);
      } else if (family_long[d] == 7) {
        real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[i]))
                   : ((flag_resid_dim == 1) ? sigma_marker_specific[i][d]
                      : sigma_y_shared[i]);
        real alpha = (P_alpha > 0) ? dot_product(X_alpha_obs[n], beta_alpha[i])
                    : alpha_skew_marker[i][d];
        target += skew_normal_lpdf(y_real[n] | eta, sig, alpha);
      } else if (family_long[d] == 8) {
        real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[i]))
                   : ((flag_resid_dim == 1) ? sigma_marker_specific[i][d]
                      : sigma_y_shared[i]);
        target += double_exponential_lpdf(y_real[n] | eta, sig);
      } else if (family_long[d] == 9) {
        real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[i]))
                   : ((flag_resid_dim == 1) ? sigma_marker_specific[i][d]
                      : sigma_y_shared[i]);
        real tau_sde = (P_tau_sde > 0) ? inv_logit(dot_product(X_tau_sde_obs[n], beta_tau_sde[i]))
                       : tau_sde_marker[i][d];
        target += skew_double_exponential_lpdf(y_real[n] | eta, sig, tau_sde);
      } else if (family_long[d] == 10) {
        real phi_beta = (P_phi_beta > 0) ? exp(dot_product(X_phi_beta_obs[n], beta_phi_beta[i]))
                        : phi_beta_marker[i][d];
        real mu = inv_logit(eta);
        real shape1 = fmax(mu * phi_beta, 1e-6);
        real shape2 = fmax((1 - mu) * phi_beta, 1e-6);
        target += beta_lpdf(y_real[n] | shape1, shape2);
      } else {
        target += ordered_logistic_lpmf(y_int[n] | eta, cutpoints_ord[i]);
      }
    }
    
    // ------------------------------------------
    // Survival Likelihood (Conditioning)
    // ------------------------------------------
    vector[15] cvm, cvk, cvm_f, cvk_f;
    
    vector[n_random_marker] vbar;
    if (n_random_marker > 0) {
      for (r in 1 : n_random_marker) {
        real acc = 0;
        for (d in 1 : n_marker_types) 
          acc += v_marker[d][r];
        vbar[r] = acc / n_marker_types;
      }
    }
    vector[n_random_marker_id] zbar;
    for (q in 1 : n_random_marker_id) {
      real acc = 0;
      for (d in 1 : n_marker_types) 
        acc += z_w[d][q];
      zbar[q] = acc / n_marker_types;
    }
    vector[n_random_marker_id] wbar_i = Li * zbar;
    
    for (j in 1 : 15) {
      cvm[j] = dot_product(mat_fixed_gk_cond[j], beta_fixed[i])
               + dot_product(mat_id_gk_cond[j], u_id);
      cvm_f[j] = dot_product(mat_fixed_gk_cond_fwd[j], beta_fixed[i])
                 + dot_product(mat_id_gk_cond_fwd[j], u_id);
      
      real mk_part = 0, mk_part_f = 0;
      if (n_random_marker > 0) {
        mk_part = dot_product(mat_marker_gk_cond[j], vbar);
        mk_part_f = dot_product(mat_marker_gk_cond_fwd[j], vbar);
      }
      cvk[j] = mk_part + dot_product(mat_marker_id_gk_cond[j], wbar_i);
      cvk_f[j] = mk_part_f
                 + dot_product(mat_marker_id_gk_cond_fwd[j], wbar_i);
    }
    
    vector[15] csm_raw = eta_fd(cvm, cvm_f, eps_finite_diff);
    vector[15] csk_raw = eta_fd(cvk, cvk_f, eps_finite_diff);
    vector[15] vcov_raw = eta_vcov_varonly_weighted_const(15, Li,
                            flag_assoc_vcov * coeff_assoc_vcov_var[i]);
    
    vector[15] cv_tot = cvm + cvk;
    vector[15] cs_tot = csm_raw + csk_raw;

    real a_cv_total = flag_assoc_cv_total * coeff_assoc_cv_total[i];
    real a_cs_total = flag_assoc_cs_total * coeff_assoc_cs_total[i];
    real a_cv_mean = flag_assoc_cv_mean * coeff_assoc_cv_mean[i];
    real a_cv_marker = flag_assoc_cv_marker * coeff_assoc_cv_marker[i];
    real a_cs_mean = flag_assoc_cs_mean * coeff_assoc_cs_mean[i];
    real a_cs_marker = flag_assoc_cs_marker * coeff_assoc_cs_marker[i];

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
      cvm,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
      coeff_cv_mean,
      spline_degree_cv_mean
    );
    vector[15] cv_marker_tf = apply_transform_vector(
      cvk,
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

    // Sum cause-specific hazards to get total hazard
    vector[15] log_h_total;
    for (j in 1 : 15) {
      vector[K_event] log_h_cause;
      for (k_ev in 1 : K_event) {
        real eta_w = 0;
        if (n_cov_hazard > 0)
          eta_w = dot_product(vec_cov_hazard, gamma_hazard[i, k_ev]);
        log_h_cause[k_ev] = log_h0_intercept[i, k_ev]
                            + dot_product(mat_basis_gk_cond[j], bs_gamma_c[i, k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    target += -cumhaz(time_condition, log_h_total, rep_vector(0.0, 15));
  }
}
generated quantities {
  #include helper/generated_quantities/dynpred_outputs.stan
}
