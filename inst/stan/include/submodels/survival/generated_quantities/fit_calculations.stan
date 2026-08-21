/**
 * @file survival/generated_quantities/fit_calculations.stan
 * @brief Calculate subject-level event likelihood, cumulative hazard and survival.
 *
 * @details
 * The declaration file introduces one output per subject.  The calculation
 * file reconstructs transformed longitudinal features separately for every
 * risk row belonging to that subject.  Summing the row contributions supports
 * ordinary right censoring, delayed entry, time-split covariates, exact events,
 * left censoring and interval censoring without changing the reported
 * subject-level form.
 *
 * For an interval-censored event in (L, R], R prepares a survival row [0, L]
 * followed by a failure-within-interval row [L, R].  The calculations below
 * therefore return log S(L) + log{1 - exp[-(H(R)-H(L))]}, which is exactly
 * log{S(L)-S(R)}.  This is a full event likelihood, not a Cox partial
 * likelihood.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  /* -------------------- Survival log-likelihood per subject */
  for (i in 1 : n_id) {
    log_lik_surv[i] = 0; // sum of all event-risk contributions belonging to subject i
    cumhaz_event[i] = 0; // cumulative hazard across all observed risk intervals for subject i

    for (e in event_start_idx[i] : event_end_idx[i]) {
    vector[n_gk] cvm_now; // Role: cvm now.
    vector[n_gk] cvm_fwd; // Role: cvm forward difference.

    for (j in 1 : n_gk) {
      cvm_now[j] = dot_product(X_gk_now[e][j], beta)
                   + dot_product(Z_id_gk_now[e][j], u_id[i]);
      cvm_fwd[j] = dot_product(X_gk_fwd[e][j], beta)
                   + dot_product(Z_id_gk_fwd[e][j], u_id[i]);
      // cvm_now/cvm_fwd: mean CV at GK node and forward shift

    }

    vector[n_gk] csm_raw = eta_fd(cvm_now, cvm_fwd, eps_fd) / tmax; // mean slope
    vector[n_gk] csk_raw; // marker slope (weighted across markers)
    int M_corr_local = num_elements(a_corr); // Role: derived dimension correlation local.
    vector[M_corr_local] corr_terms_raw = eta_chol_corr(L_i[i]); // raw off-diagonal K_i terms
    int M_vcov_local = num_elements(a_vcov); // Role: derived dimension covariance local.
    matrix[Q_idm, Q_idm] L_i_eff_assoc = L_i[i]; // subject-specific covariance factor on the original-time basis
    vector[M_vcov_local] vcov_terms_raw = eta_vcov_weighted_const(L_i_eff_assoc, (M_vcov_local == Q_idm)); // raw K_i off-diagonals plus effective SD_i entries

    int D_mkrs = size(v_marker); // Role: D mkrs.
    vector[n_gk] cv_tot_tf; // Role: current value total transformation.
    vector[n_gk] cv_mean_tf = apply_transform_vector(
      cvm_now,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      functional_iota_intercept_idx_cv_mean,
      functional_iota_slope_idx_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
      coeff_cv_mean_eff,
      spline_degree_cv_mean,
      iota_intercept_cv_mean,
      iota_slope_cv_mean
    );
    vector[n_gk] cv_marker_tf; // Role: current value marker transformation.
    vector[n_gk] cs_tot_tf; // Role: current slope total transformation.
    vector[n_gk] cs_marker_tf; // Role: current slope marker transformation.
    for (j in 1 : n_gk) {
      real acc_cv_tot_tf = 0; // Role: acc current value total transformation.
      real acc_cv_marker_tf = 0; // Role: acc current value marker transformation.
      real acc_csk_raw = 0; // Role: acc csk unscaled.
      real acc_cs_tot_tf = 0; // Role: acc current slope total transformation.
      real acc_cs_marker_tf = 0; // Role: acc current slope marker transformation.
      for (d in 1 : D_mkrs) {
        real mk_part_now_d = 0; // Role: mk part now d.
        real mk_part_fwd_d = 0; // Role: mk part forward difference d.
        if (R_mk > 0) {
          mk_part_now_d = dot_product(Z_mk_gk_now[e][j], v_marker[d]);
          mk_part_fwd_d = dot_product(Z_mk_gk_fwd[e][j], v_marker[d]);
        }
        real cvk_now_d = mk_part_now_d + dot_product(Z_idm_gk_now[e][j], w_idm[i, d]); // Role: cvk now d.
        real cvk_fwd_d = mk_part_fwd_d + dot_product(Z_idm_gk_fwd[e][j], w_idm[i, d]); // Role: cvk forward difference d.
        real cv_tot_now_d = cvm_now[j] + cvk_now_d; // Role: current value total now d.
        acc_cv_tot_tf += marker_weights_eff_cv_total[d]
          * apply_transform_scalar(cv_tot_now_d, tf_mode_cv_tot,
                                   functional_ops_cv, functional_iota_intercept_idx_cv, functional_iota_slope_idx_cv, const_data_cv,
                                   knots_cv, coeff_cv_eff, spline_degree_cv,
                                   iota_intercept_cv, iota_slope_cv);
        acc_cv_marker_tf += marker_weights_eff_cv_marker[d]
          * apply_transform_scalar(cvk_now_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, functional_iota_intercept_idx_cv_marker, functional_iota_slope_idx_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker_eff, spline_degree_cv_marker,
                                   iota_intercept_cv_marker, iota_slope_cv_marker);
        {
          real csk_raw_d = ((cvk_fwd_d - cvk_now_d) / eps_fd) / tmax; // Role: csk unscaled d.
          real cs_tot_raw_d = csm_raw[j] + csk_raw_d; // Role: current slope total unscaled d.

          acc_csk_raw += marker_weights_eff_cs_marker[d] * csk_raw_d;
          acc_cs_tot_tf += marker_weights_eff_cs_total[d]
            * apply_transform_scalar(cs_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, functional_iota_intercept_idx_cs, functional_iota_slope_idx_cs, const_data_cs,
                                     knots_cs, coeff_cs_eff, spline_degree_cs,
                                     iota_intercept_cs, iota_slope_cs);
          acc_cs_marker_tf += marker_weights_eff_cs_marker[d]
            * apply_transform_scalar(csk_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, functional_iota_intercept_idx_cs_marker, functional_iota_slope_idx_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker_eff, spline_degree_cs_marker,
                                     iota_intercept_cs_marker, iota_slope_cs_marker);
        }
      }
      cv_tot_tf[j] = acc_cv_tot_tf / D_mkrs;
      cv_marker_tf[j] = acc_cv_marker_tf / D_mkrs;
      csk_raw[j] = acc_csk_raw / D_mkrs;
      cs_tot_tf[j] = acc_cs_tot_tf / D_mkrs;
      cs_marker_tf[j] = acc_cs_marker_tf / D_mkrs;
    }

    vector[n_gk] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      functional_iota_intercept_idx_cs_mean,
      functional_iota_slope_idx_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
      coeff_cs_mean_eff,
      spline_degree_cs_mean,
      iota_intercept_cs_mean,
      iota_slope_cs_mean
    );
    vector[M_corr_local] corr_terms_tf = apply_transform_vector_by_component(
      corr_terms_raw,
      tf_mode_corr,
      functional_ops_corr,
      functional_iota_intercept_idx_corr,
      functional_iota_slope_idx_corr,
      const_data_corr,
      knots_corr,
      coeff_corr_eff,
      spline_degree_corr,
      iota_intercept_corr,
      iota_slope_corr
    );
    vector[M_corr_local] corr_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_corr_local),
      tf_mode_corr,
      functional_ops_corr,
      functional_iota_intercept_idx_corr,
      functional_iota_slope_idx_corr,
      const_data_corr,
      knots_corr,
      coeff_corr_eff,
      spline_degree_corr,
      iota_intercept_corr,
      iota_slope_corr
    );
    corr_terms_tf = corr_terms_tf - corr_terms_ref;
    real corr_assoc_scalar = dot_product(a_corr, corr_terms_tf); // Role: correlation assoc scalar.
    vector[n_gk] corr_assoc = rep_vector(corr_assoc_scalar, n_gk); // Role: correlation assoc.
    vector[M_vcov_local] vcov_terms_tf = apply_transform_vector_by_component(
      vcov_terms_raw,
      tf_mode_vcov,
      functional_ops_vcov,
      functional_iota_intercept_idx_vcov,
      functional_iota_slope_idx_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov_eff,
      spline_degree_vcov,
      iota_intercept_vcov,
      iota_slope_vcov
    );
    vector[M_vcov_local] vcov_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_vcov_local),
      tf_mode_vcov,
      functional_ops_vcov,
      functional_iota_intercept_idx_vcov,
      functional_iota_slope_idx_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov_eff,
      spline_degree_vcov,
      iota_intercept_vcov,
      iota_slope_vcov
    );
    vcov_terms_tf = vcov_terms_tf - vcov_terms_ref;
    real vcov_assoc_scalar = dot_product(a_vcov, vcov_terms_tf); // Role: covariance assoc scalar.
    vector[n_gk] vcov_assoc = rep_vector(vcov_assoc_scalar, n_gk); // Role: covariance assoc.

    vector[n_gk] eta_assoc_nodes = a_cv_total * cv_tot_tf
                                 + a_cv_mean * cv_mean_tf
                                 + a_cv_marker * cv_marker_tf
                                 + a_cs_total * cs_tot_tf
                                 + a_cs_mean * cs_mean_tf
                                 + a_cs_marker * cs_marker_tf
                                 + corr_assoc
                                 + vcov_assoc;
    // eta_assoc_nodes: association predictor across GK nodes

    // Event contribution
    real event_term = 0; // log-hazard at event time (if event)
    if (event_censor_type[e] == 1) {
      int k_ev = event_type[e]; // cause of the exactly observed event at this row's upper limit
      real eta_w_ev = 0; // hazard covariate contribution
      if (p_w > 0)
        eta_w_ev = dot_product(to_vector(W[e]'), gamma_w[k_ev]);
      real log_h0_S = dot_product(Bs_event_c[e], bs_gamma_c[k_ev]); // Role: log h0 S.
      // log_h0_S: baseline log-hazard at event time

      real cvm_S = dot_product(X_event_now[e], beta)
                   + dot_product(Z_id_event_now[e], u_id[i]);
      real cvm_S_fwd = dot_product(X_event_fwd[e], beta)
                       + dot_product(Z_id_event_fwd[e], u_id[i]);
      // cvm_S/cvm_S_fwd: mean CV at event time and forward shift
      real mk_part_S = 0; // Role: mk part S.
      real mk_part_S_fwd = 0; // Role: mk part S forward difference.
      // mk_part_S*: marker-only CV at event time
      if (R_mk > 0) {
        mk_part_S = dot_product(Z_mk_event_now[e], vbar);
        mk_part_S_fwd = dot_product(Z_mk_event_fwd[e], vbar);
      }
      real cvk_S = mk_part_S + dot_product(Z_idm_event_now[e], wbar_i[i]); // Role: cvk S.
      real cvk_S_fwd = mk_part_S_fwd
               + dot_product(Z_idm_event_fwd[e], wbar_i[i]);
      // cvk_S/cvk_S_fwd: marker-id CV at event time and forward shift

      real csm_S_raw = ((cvm_S_fwd - cvm_S) / eps_fd) / tmax; // mean slope at event
      real csk_S_raw = 0; // marker slope at event (weighted across markers)

      real corr_S_assoc = corr_assoc_scalar; // same transformed-and-weighted constant at event time
      real vcov_S_assoc = vcov_assoc_scalar; // same transformed-and-weighted covariance contribution at event time

      real cv_S_tot_tf = 0; // Role: current value S total transformation.
      real cv_S_mean_tf = apply_transform_scalar(
        cvm_S,
        tf_mode_cv_mean,
        functional_ops_cv_mean,
        functional_iota_intercept_idx_cv_mean,
        functional_iota_slope_idx_cv_mean,
        const_data_cv_mean,
        knots_cv_mean,
        coeff_cv_mean_eff,
        spline_degree_cv_mean,
        iota_intercept_cv_mean,
        iota_slope_cv_mean
      );
      real cv_S_marker_tf = 0; // Role: current value S marker transformation.
      real cs_S_tot_tf = 0; // Role: current slope S total transformation.
      real cs_S_marker_tf = 0; // Role: current slope S marker transformation.
      for (d in 1 : D_mkrs) {
        real mk_part_S_d = 0; // Role: mk part S d.
        real mk_part_S_fwd_d = 0; // Role: mk part S forward difference d.
        if (R_mk > 0) {
          mk_part_S_d = dot_product(Z_mk_event_now[e], v_marker[d]);
          mk_part_S_fwd_d = dot_product(Z_mk_event_fwd[e], v_marker[d]);
        }
        real cvk_S_d = mk_part_S_d + dot_product(Z_idm_event_now[e], w_idm[i, d]); // Role: cvk S d.
        real cvk_S_fwd_d = mk_part_S_fwd_d + dot_product(Z_idm_event_fwd[e], w_idm[i, d]); // Role: cvk S forward difference d.
        real cv_S_tot_d = cvm_S + cvk_S_d; // Role: current value S total d.

        cv_S_tot_tf += marker_weights_eff_cv_total[d]
          * apply_transform_scalar(cv_S_tot_d, tf_mode_cv_tot,
                                   functional_ops_cv, functional_iota_intercept_idx_cv, functional_iota_slope_idx_cv, const_data_cv,
                                   knots_cv, coeff_cv_eff, spline_degree_cv,
                                   iota_intercept_cv, iota_slope_cv);
        cv_S_marker_tf += marker_weights_eff_cv_marker[d]
          * apply_transform_scalar(cvk_S_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, functional_iota_intercept_idx_cv_marker, functional_iota_slope_idx_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker_eff, spline_degree_cv_marker,
                                   iota_intercept_cv_marker, iota_slope_cv_marker);
        {
          real csk_S_raw_d = ((cvk_S_fwd_d - cvk_S_d) / eps_fd) / tmax; // Role: csk S unscaled d.
          real cs_S_tot_raw_d = csm_S_raw + csk_S_raw_d; // Role: current slope S total unscaled d.

          csk_S_raw += marker_weights_eff_cs_marker[d] * csk_S_raw_d;
          cs_S_tot_tf += marker_weights_eff_cs_total[d]
            * apply_transform_scalar(cs_S_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, functional_iota_intercept_idx_cs, functional_iota_slope_idx_cs, const_data_cs,
                                     knots_cs, coeff_cs_eff, spline_degree_cs,
                                     iota_intercept_cs, iota_slope_cs);
          cs_S_marker_tf += marker_weights_eff_cs_marker[d]
            * apply_transform_scalar(csk_S_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, functional_iota_intercept_idx_cs_marker, functional_iota_slope_idx_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker_eff, spline_degree_cs_marker,
                                     iota_intercept_cs_marker, iota_slope_cs_marker);
        }
      }
      cv_S_tot_tf /= D_mkrs;
      cv_S_marker_tf /= D_mkrs;
      csk_S_raw /= D_mkrs;
      cs_S_tot_tf /= D_mkrs;
      cs_S_marker_tf /= D_mkrs;

      real cs_S_mean_tf = apply_transform_scalar(
        csm_S_raw,
        tf_mode_cs_mean,
        functional_ops_cs_mean,
        functional_iota_intercept_idx_cs_mean,
        functional_iota_slope_idx_cs_mean,
        const_data_cs_mean,
        knots_cs_mean,
        coeff_cs_mean_eff,
        spline_degree_cs_mean,
        iota_intercept_cs_mean,
        iota_slope_cs_mean
      );
      real eta_assoc_S = a_cv_total * cv_S_tot_tf
                         + a_cv_mean * cv_S_mean_tf
                         + a_cv_marker * cv_S_marker_tf
                         + a_cs_total * cs_S_tot_tf
                         + a_cs_mean * cs_S_mean_tf
                         + a_cs_marker * cs_S_marker_tf
                         + corr_S_assoc
                         + vcov_S_assoc;
      // eta_assoc_S: association predictor at event time

      event_term = log_h0_S + eta_w_ev + eta_assoc_S;
    }

    vector[n_gk] log_h_total; // total log-hazard at quadrature nodes
    for (j in 1 : n_gk) {
      vector[K_event] log_h_cause; // cause-specific log hazards
      for (k_ev in 1 : K_event) {
        real eta_w = 0; // Role: linear predictor w.
        if (p_w > 0)
          eta_w = dot_product(to_vector(W_gk[e][j]'), gamma_w[k_ev]);
        log_h_cause[k_ev] = dot_product(Bs_gk_c[e][j], bs_gamma_c[k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    {
      real interval_width = S_event[e] - S_entry[e]; // scaled duration represented by event-risk row e
      real cumulative_hazard_interval = cumhaz(interval_width, n_gk, log_h_total, rep_vector(0.0, n_gk)); // integrated all-cause hazard over row e
      cumhaz_event[i] += cumulative_hazard_interval;

      if (event_censor_type[e] == 0) {
        log_lik_surv[i] += -cumulative_hazard_interval; // known survival throughout this interval
      } else if (event_censor_type[e] == 1) {
        log_lik_surv[i] += event_term - cumulative_hazard_interval; // exact failure density at the upper limit
      } else {
        log_lik_surv[i] += log_failure_within_interval(cumulative_hazard_interval); // failure somewhere within the interval
      }
    }
    } // finish the event-risk rows belonging to subject i

    surv_prob_event[i] = exp(-cumhaz_event[i]); // survival through all represented risk intervals for subject i
  }
