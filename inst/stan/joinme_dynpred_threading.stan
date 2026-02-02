/**
 * @file joinme_dynpred_threading.stan
 * @brief Threaded dynamic prediction for joint longitudinal + survival models.
 *
 * ## Threading strategy
 * This variant parallelizes *draw-wise* random-effects sampling and the
 * conditioning likelihood using `reduce_sum`. Each task operates on a slice
 * of posterior draws, which is typically cheaper to balance than subject-wise
 * splits during prediction. The `grainsize` parameter controls how many draws
 * each task receives (smaller = more tasks, better load balance; larger = less
 * overhead and sometimes better cache locality).
 *
 * ## Model parity
 * The prediction logic matches the non-threaded version exactly:
 * - Same random-effects reconstruction and covariance regression.
 * - Same association term construction and transform pipeline.
 * - Same fitted values, predictive draws, and survival probabilities.
 *
 * Note: No RNG is used inside the reduce_sum partial; all randomness comes from
 * posterior draws already supplied as data. This keeps results deterministic
 * across thread counts.
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/functions/functional_transform.stanfunctions
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions

  /**
   * @brief Per-draw log-likelihood slice for reduce_sum.
   *
   * Each call evaluates a contiguous slice of posterior draws (draw_ids_slice)
   * and returns the summed log-likelihood contributions for those draws.
   * The logic mirrors the non-threaded model block, but accumulates into `lp`
   * instead of `target` because this is a pure function.
   */
  real partial_draw(
      array[] int draw_ids_slice,
      int start,
      int end,
      int n_obs_long,
      array[] int idx_marker_obs,
      int n_marker_types,
      vector y_real,
      array[] int y_int,
      array[] int trials_obs,
      int n_fixed_effects,
      int n_random_id,
      int n_random_marker,
      int n_random_marker_id,
      matrix mat_fixed_obs,
      matrix mat_id_obs,
      matrix mat_marker_obs,
      matrix mat_marker_id_obs,
      int P_sigma,
      matrix X_sigma_obs,
      int P_nu,
      matrix X_nu_obs,
      int P_phi,
      matrix X_phi_obs,
      int P_alpha,
      matrix X_alpha_obs,
      int P_phi_beta,
      matrix X_phi_beta_obs,
      int P_tau_sde,
      matrix X_tau_sde_obs,
      int n_cov_vcov,
      vector vec_cov_vcov,
      int n_cov_hazard,
      vector vec_cov_hazard,
      int n_basehaz_basis,
      real time_condition,
      matrix mat_basis_gk_cond,
      matrix mat_fixed_gk_cond,
      matrix mat_id_gk_cond,
      matrix mat_marker_gk_cond,
      matrix mat_marker_id_gk_cond,
      matrix mat_fixed_gk_cond_fwd,
      matrix mat_id_gk_cond_fwd,
      matrix mat_marker_gk_cond_fwd,
      matrix mat_marker_id_gk_cond_fwd,
      real eps_finite_diff,
      array[] int family_long,
      int flag_resid_dim,
      int num_unique_cov_entries,
      array[] int idx_row_cov,
      array[] int idx_col_cov,
      array[] vector beta_fixed,
      array[] vector beta_sigma,
      array[] vector beta_nu,
      array[] vector beta_phi,
      array[] vector beta_alpha,
      array[] vector beta_phi_beta,
      array[] vector beta_tau_sde,
      array[] vector tau_id,
      array[] matrix Lcorr_id,
      array[] vector tau_marker,
      array[] matrix Lcorr_marker,
      array[] vector tau_marker_id,
      array[] matrix Lcorr_marker_id,
      array[] matrix B_cross,
      array[] vector alpha_vcov_reg,
      array[] vector beta_vcov_reg_flat,
      array[] real tau_vcov_reg,
      array[] vector lambda_vcov_reg,
      int K_event,
      array[,] real log_h0_intercept,
      array[,] vector bs_gamma_c,
      array[,] vector gamma_hazard,
      array[] real sigma_y_shared,
      array[] vector sigma_marker_specific,
      array[] vector nu_marker,
      array[] vector phi_nb_marker,
      array[] vector alpha_skew_marker,
      array[] vector phi_beta_marker,
      array[] vector tau_sde_marker,
      array[] real coeff_assoc_cv_total,
      array[] real coeff_assoc_cs_total,
      array[] real coeff_assoc_cv_mean,
      array[] real coeff_assoc_cs_mean,
      array[] real coeff_assoc_cv_marker,
      array[] real coeff_assoc_cs_marker,
      array[] vector coeff_assoc_vcov_var,
      int K_ord,
      array[] vector cutpoints_ord,
      int flag_assoc_cv_total,
      int flag_assoc_cs_total,
      int flag_assoc_cv_mean,
      int flag_assoc_cv_marker,
      int flag_assoc_cs_mean,
      int flag_assoc_cs_marker,
      int flag_assoc_vcov,
      int tf_mode_cv_tot,
      int tf_mode_cs_tot,
      int tf_mode_vcov,
      int tf_mode_cv_mean,
      int tf_mode_cv_marker,
      int tf_mode_cs_mean,
      int tf_mode_cs_marker,
      array[] int functional_ops_cv,
      vector const_data_cv,
      vector knots_cv,
      vector coeff_cv,
      int spline_degree_cv,
      array[] int functional_ops_cs,
      vector const_data_cs,
      vector knots_cs,
      vector coeff_cs,
      int spline_degree_cs,
      array[] int functional_ops_vcov,
      vector const_data_vcov,
      vector knots_vcov,
      vector coeff_vcov,
      int spline_degree_vcov,
      array[] int functional_ops_cv_mean,
      vector const_data_cv_mean,
      vector knots_cv_mean,
      vector coeff_cv_mean,
      int spline_degree_cv_mean,
      array[] int functional_ops_cv_marker,
      vector const_data_cv_marker,
      vector knots_cv_marker,
      vector coeff_cv_marker,
      int spline_degree_cv_marker,
      array[] int functional_ops_cs_mean,
      vector const_data_cs_mean,
      vector knots_cs_mean,
      vector coeff_cs_mean,
      int spline_degree_cs_mean,
      array[] int functional_ops_cs_marker,
      vector const_data_cs_marker,
      vector knots_cs_marker,
      vector coeff_cs_marker,
      int spline_degree_cs_marker,
      int flag_indep_id_re,
      int flag_indep_marker_re,
      int flag_indep_marker_byid_latent_re,
      int flag_allow_marker_crosscorr,
      array[] vector z_u,
      array[,] vector z_v,
      array[,] vector z_w_lat,
      vector z_L
  ) {
    real lp = 0;

    for (i in start:end) {
      int k = draw_ids_slice[i];

      lp += std_normal_lpdf(z_u[k]);
      if (n_random_marker > 0) {
        for (d in 1 : n_marker_types)
          lp += std_normal_lpdf(z_v[k, d]);
      }
      for (d in 1 : n_marker_types)
        lp += std_normal_lpdf(z_w_lat[k, d]);
      lp += std_normal_lpdf(z_L[k]);

      matrix[n_random_id, n_random_id] L_u;
      if (flag_indep_id_re == 1)
        L_u = diag_matrix(tau_id[k]);
      else
        L_u = diag_pre_multiply(tau_id[k], Lcorr_id[k]);
      vector[n_random_id] u_id = L_u * z_u[k];

      matrix[n_random_marker, n_random_marker] L_v;
      if (n_random_marker > 0) {
        if (flag_indep_marker_re == 1)
          L_v = diag_matrix(tau_marker[k]);
        else
          L_v = diag_pre_multiply(tau_marker[k], Lcorr_marker[k]);
      }
      array[n_marker_types] vector[n_random_marker] v_marker;
      if (n_random_marker > 0) {
        for (d in 1 : n_marker_types)
          v_marker[d] = L_v * z_v[k, d];
      }

      matrix[n_random_marker_id, n_random_marker_id] L_w;
      if (flag_indep_marker_byid_latent_re == 1)
        L_w = diag_matrix(tau_marker_id[k]);
      else
        L_w = diag_pre_multiply(tau_marker_id[k], Lcorr_marker_id[k]);
      array[n_marker_types] vector[n_random_marker_id] z_w;
      for (d in 1 : n_marker_types) {
        vector[n_random_marker_id] cross = rep_vector(0.0, n_random_marker_id);
        if (flag_allow_marker_crosscorr == 1 && n_random_marker > 0)
          cross = B_cross[k] * v_marker[d];
        z_w[d] = cross + L_w * z_w_lat[k, d];
      }

      real u_Lk = tau_vcov_reg[k] * z_L[k];
      matrix[n_random_marker_id, n_random_marker_id] Li = rep_matrix(0.0,
                                                                     n_random_marker_id,
                                                                     n_random_marker_id);
      for (m in 1 : num_unique_cov_entries) {
        int r_ = idx_row_cov[m];
        int c_ = idx_col_cov[m];
        vector[n_cov_vcov] bL_m = beta_vcov_reg_flat[k][((m - 1) * n_cov_vcov
                                                         + 1) : (m * n_cov_vcov)];
        real lp_i = alpha_vcov_reg[k][m] + dot_product(bL_m, vec_cov_vcov)
                    + lambda_vcov_reg[k][m] * u_Lk;
        Li[r_, c_] = (r_ == c_) ? log1p_exp(lp_i) : lp_i;
      }

      array[n_marker_types] vector[n_random_marker_id] w_idscaled;
      for (d in 1 : n_marker_types)
        w_idscaled[d] = Li * z_w[d];

      for (n in 1 : n_obs_long) {
        int d = idx_marker_obs[n];
        real eta = dot_product(mat_fixed_obs[n], beta_fixed[k])
                   + dot_product(mat_id_obs[n], u_id)
                   + ((n_random_marker > 0)
                      ? dot_product(mat_marker_obs[n], v_marker[d]) : 0.0)
                   + dot_product(mat_marker_id_obs[n], w_idscaled[d]);

        if (family_long[d] == 1) {
          real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[k]))
                     : ((flag_resid_dim == 1) ? sigma_marker_specific[k][d]
                        : sigma_y_shared[k]);
          lp += normal_lpdf(y_real[n] | eta, sig);
        } else if (family_long[d] == 2) {
          real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[k]))
                     : ((flag_resid_dim == 1) ? sigma_marker_specific[k][d]
                        : sigma_y_shared[k]);
          real nu = (P_nu > 0) ? (2 + exp(dot_product(X_nu_obs[n], beta_nu[k])))
                    : nu_marker[k][d];
          lp += student_t_lpdf(y_real[n] | nu, eta, sig);
        } else if (family_long[d] == 3) {
          lp += bernoulli_logit_lpmf(y_int[n] | eta);
        } else if (family_long[d] == 4) {
          lp += binomial_logit_lpmf(y_int[n] | trials_obs[n], eta);
        } else if (family_long[d] == 5) {
          lp += poisson_log_lpmf(y_int[n] | eta);
        } else if (family_long[d] == 6) {
          real phi = (P_phi > 0) ? exp(dot_product(X_phi_obs[n], beta_phi[k]))
                     : phi_nb_marker[k][d];
          lp += neg_binomial_2_log_lpmf(y_int[n] | eta, phi);
        } else if (family_long[d] == 7) {
          real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[k]))
                     : ((flag_resid_dim == 1) ? sigma_marker_specific[k][d]
                        : sigma_y_shared[k]);
          real alpha = (P_alpha > 0) ? dot_product(X_alpha_obs[n], beta_alpha[k])
                      : alpha_skew_marker[k][d];
          lp += skew_normal_lpdf(y_real[n] | eta, sig, alpha);
        } else if (family_long[d] == 8) {
          real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[k]))
                     : ((flag_resid_dim == 1) ? sigma_marker_specific[k][d]
                        : sigma_y_shared[k]);
          lp += double_exponential_lpdf(y_real[n] | eta, sig);
        } else if (family_long[d] == 9) {
          real sig = (P_sigma > 0) ? exp(dot_product(X_sigma_obs[n], beta_sigma[k]))
                     : ((flag_resid_dim == 1) ? sigma_marker_specific[k][d]
                        : sigma_y_shared[k]);
          real tau_sde = (P_tau_sde > 0) ? inv_logit(dot_product(X_tau_sde_obs[n], beta_tau_sde[k]))
                         : tau_sde_marker[k][d];
          lp += skew_double_exponential_lpdf(y_real[n] | eta, sig, tau_sde);
        } else if (family_long[d] == 10) {
          real phi_beta = (P_phi_beta > 0) ? exp(dot_product(X_phi_beta_obs[n], beta_phi_beta[k]))
                          : phi_beta_marker[k][d];
          real mu = inv_logit(eta);
          real shape1 = fmax(mu * phi_beta, 1e-6);
          real shape2 = fmax((1 - mu) * phi_beta, 1e-6);
          lp += beta_lpdf(y_real[n] | shape1, shape2);
        } else {
          lp += ordered_logistic_lpmf(y_int[n] | eta, cutpoints_ord[k]);
        }
      }

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
        cvm[j] = dot_product(mat_fixed_gk_cond[j], beta_fixed[k])
                 + dot_product(mat_id_gk_cond[j], u_id);
        cvm_f[j] = dot_product(mat_fixed_gk_cond_fwd[j], beta_fixed[k])
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
                              flag_assoc_vcov * coeff_assoc_vcov_var[k]);

      vector[15] cv_tot = cvm + cvk;
      vector[15] cs_tot = csm_raw + csk_raw;

      real a_cv_total = flag_assoc_cv_total * coeff_assoc_cv_total[k];
      real a_cs_total = flag_assoc_cs_total * coeff_assoc_cs_total[k];
      real a_cv_mean = flag_assoc_cv_mean * coeff_assoc_cv_mean[k];
      real a_cv_marker = flag_assoc_cv_marker * coeff_assoc_cv_marker[k];
      real a_cs_mean = flag_assoc_cs_mean * coeff_assoc_cs_mean[k];
      real a_cs_marker = flag_assoc_cs_marker * coeff_assoc_cs_marker[k];

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

      vector[15] log_h_total;
      for (j in 1 : 15) {
        vector[K_event] log_h_cause;
        for (k_ev in 1 : K_event) {
          real eta_w = 0;
          if (n_cov_hazard > 0)
            eta_w = dot_product(vec_cov_hazard, gamma_hazard[k, k_ev]);
          log_h_cause[k_ev] = log_h0_intercept[k, k_ev]
                              + dot_product(mat_basis_gk_cond[j], bs_gamma_c[k, k_ev])
                              + eta_w
                              + eta_assoc_nodes[j];
        }
        log_h_total[j] = log_sum_exp(log_h_cause);
      }

      lp += -cumhaz(time_condition, log_h_total, rep_vector(0.0, 15));
    }

    return lp;
  }
}

data {
  #include helper/data/joinme_dynpred_common.stan
  int<lower=1> grainsize;

  int flag_indep_id_re;
  int flag_indep_marker_re;
  int flag_indep_marker_byid_latent_re;
  int flag_indep_idmarker_cov;
  int flag_allow_marker_crosscorr;

  array[num_unique_cov_entries] int idx_row_cov;
  array[num_unique_cov_entries] int idx_col_cov;
}

transformed data {
  array[n_draws] int draw_ids;
  for (i in 1 : n_draws) draw_ids[i] = i;
}

parameters {
  #include helper/parameters/joinme_dynpred_common.stan
}

model {
  // -------------------- Threaded draw-wise conditioning likelihood
  // Each task handles a slice of draws; this is typically more stable than
  // subject-wise splits for prediction workloads.
  target += reduce_sum(
    partial_draw,
    draw_ids,
    grainsize,
    n_obs_long,
    idx_marker_obs,
    n_marker_types,
    y_real,
    y_int,
    trials_obs,
    n_fixed_effects,
    n_random_id,
    n_random_marker,
    n_random_marker_id,
    mat_fixed_obs,
    mat_id_obs,
    mat_marker_obs,
    mat_marker_id_obs,
    P_sigma,
    X_sigma_obs,
    P_nu,
    X_nu_obs,
    P_phi,
    X_phi_obs,
    P_alpha,
    X_alpha_obs,
    P_phi_beta,
    X_phi_beta_obs,
    P_tau_sde,
    X_tau_sde_obs,
    n_cov_vcov,
    vec_cov_vcov,
    n_cov_hazard,
    vec_cov_hazard,
    n_basehaz_basis,
    time_condition,
    mat_basis_gk_cond,
    mat_fixed_gk_cond,
    mat_id_gk_cond,
    mat_marker_gk_cond,
    mat_marker_id_gk_cond,
    mat_fixed_gk_cond_fwd,
    mat_id_gk_cond_fwd,
    mat_marker_gk_cond_fwd,
    mat_marker_id_gk_cond_fwd,
    eps_finite_diff,
    family_long,
    flag_resid_dim,
    num_unique_cov_entries,
    idx_row_cov,
    idx_col_cov,
    beta_fixed,
    beta_sigma,
    beta_nu,
    beta_phi,
    beta_alpha,
    beta_phi_beta,
    beta_tau_sde,
    tau_id,
    Lcorr_id,
    tau_marker,
    Lcorr_marker,
    tau_marker_id,
    Lcorr_marker_id,
    B_cross,
    alpha_vcov_reg,
    beta_vcov_reg_flat,
    tau_vcov_reg,
    lambda_vcov_reg,
    K_event,
    log_h0_intercept,
    bs_gamma_c,
    gamma_hazard,
    sigma_y_shared,
    sigma_marker_specific,
    nu_marker,
    phi_nb_marker,
    alpha_skew_marker,
    phi_beta_marker,
    tau_sde_marker,
    coeff_assoc_cv_total,
    coeff_assoc_cs_total,
    coeff_assoc_cv_mean,
    coeff_assoc_cs_mean,
    coeff_assoc_cv_marker,
    coeff_assoc_cs_marker,
    coeff_assoc_vcov_var,
    K_ord,
    cutpoints_ord,
    flag_assoc_cv_total,
    flag_assoc_cs_total,
    flag_assoc_cv_mean,
    flag_assoc_cv_marker,
    flag_assoc_cs_mean,
    flag_assoc_cs_marker,
    flag_assoc_vcov,
    tf_mode_cv_tot,
    tf_mode_cs_tot,
    tf_mode_vcov,
    tf_mode_cv_mean,
    tf_mode_cv_marker,
    tf_mode_cs_mean,
    tf_mode_cs_marker,
    functional_ops_cv,
    const_data_cv,
    knots_cv,
    coeff_cv,
    spline_degree_cv,
    functional_ops_cs,
    const_data_cs,
    knots_cs,
    coeff_cs,
    spline_degree_cs,
    functional_ops_vcov,
    const_data_vcov,
    knots_vcov,
    coeff_vcov,
    spline_degree_vcov,
    functional_ops_cv_mean,
    const_data_cv_mean,
    knots_cv_mean,
    coeff_cv_mean,
    spline_degree_cv_mean,
    functional_ops_cv_marker,
    const_data_cv_marker,
    knots_cv_marker,
    coeff_cv_marker,
    spline_degree_cv_marker,
    functional_ops_cs_mean,
    const_data_cs_mean,
    knots_cs_mean,
    coeff_cs_mean,
    spline_degree_cs_mean,
    functional_ops_cs_marker,
    const_data_cs_marker,
    knots_cs_marker,
    coeff_cs_marker,
    spline_degree_cs_marker,
    flag_indep_id_re,
    flag_indep_marker_re,
    flag_indep_marker_byid_latent_re,
    flag_allow_marker_crosscorr,
    z_u,
    z_v,
    z_w_lat,
    z_L
  );
}

generated quantities {
  #include helper/generated_quantities/dynpred_outputs.stan
}
