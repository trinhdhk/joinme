/**
 * @file functional/model/fit.stan
 * @brief Assign shape and smoothness priors to association transformations.
 *
 * @details
 * Piecewise-linear increments are exchangeable on their simplexes.  I-spline increments follow the association slope family's raw law.  Optional second-difference penalties discourage unnecessary curvature without altering monotonicity.
 */
  /**
   * @brief Exchangeable Dirichlet priors for piecewise-linear increments.
   *
   * @details A vector of ones gives a uniform prior over the simplex, so no
   * interval is preferred before seeing the longitudinal and event data. The
   * optional second-difference penalty below can add smoothness without
   * changing the ordering constraint. Existing I-splines retain latent
   * increment priors chosen by the association-slope family below.
   */
  if ((tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) && n_free_spline_cv > 0)
    pwlin_simplex_cv ~ dirichlet(rep_vector(1, n_free_spline_cv));
  if ((tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) && n_free_spline_cs > 0)
    pwlin_simplex_cs ~ dirichlet(rep_vector(1, n_free_spline_cs));
  if ((tf_mode_corr == 3 || tf_mode_corr == 7) && n_free_spline_corr > 0)
    for (m in 1:M_corr)
      pwlin_simplex_corr[m] ~ dirichlet(rep_vector(1, n_free_spline_corr));
  if ((tf_mode_vcov == 3 || tf_mode_vcov == 7) && n_free_spline_vcov > 0)
    for (m in 1:M_vcov)
      pwlin_simplex_vcov[m] ~ dirichlet(rep_vector(1, n_free_spline_vcov));
  if ((tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) && n_free_spline_cv_mean > 0)
    pwlin_simplex_cv_mean ~ dirichlet(rep_vector(1, n_free_spline_cv_mean));
  if ((tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) && n_free_spline_cv_marker > 0)
    pwlin_simplex_cv_marker ~ dirichlet(rep_vector(1, n_free_spline_cv_marker));
  if ((tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) && n_free_spline_cs_mean > 0)
    pwlin_simplex_cs_mean ~ dirichlet(rep_vector(1, n_free_spline_cs_mean));
  if ((tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) && n_free_spline_cs_marker > 0)
    pwlin_simplex_cs_marker ~ dirichlet(rep_vector(1, n_free_spline_cs_marker));

  /*
   * Association-shape increments use the association-slope family's standard
   * raw law.
   * A regularised horseshoe has a standard-Normal raw law, so the simplex or
   * softmax shape coordinates remain separate from the horseshoe coefficient
   * scales and do not duplicate their global--local hierarchy.
   */
  if (prior_regression_family[prior_start_alpha] == 3) {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ double_exponential(0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ double_exponential(0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ double_exponential(0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ double_exponential(0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ double_exponential(0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ double_exponential(0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ double_exponential(0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ double_exponential(0, 1);

  } else if (
    prior_regression_family[prior_start_alpha] == 2 ||
    prior_regression_family[prior_start_alpha] == 4
  ) {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ std_normal();
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ std_normal();
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ std_normal();
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ std_normal();
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ std_normal();
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ std_normal();
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ std_normal();
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ std_normal();

  } else {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ student_t(prior_regression_df[prior_start_alpha], 0, 1);

  }

  

  // - first differences control monotone increase,
  // - second differences control wiggliness,
  // - lambda says how strongly we discourage bends.
  if (estimate_spline_cv == 1 && n_coeff_cv > 2 && lambda_spline_cv > 0) {
    for (j in 3:n_coeff_cv)
      target += -0.5 * lambda_spline_cv * square(coeff_cv_eff[j] - 2 * coeff_cv_eff[j - 1] + coeff_cv_eff[j - 2]);
  }
  if (estimate_spline_cs == 1 && n_coeff_cs > 2 && lambda_spline_cs > 0) {
    for (j in 3:n_coeff_cs)
      target += -0.5 * lambda_spline_cs * square(coeff_cs_eff[j] - 2 * coeff_cs_eff[j - 1] + coeff_cs_eff[j - 2]);
  }
  if (estimate_spline_corr == 1 && n_coeff_corr > 2 && lambda_spline_corr > 0) {
    for (m in 1:M_corr)
      for (j in 3:n_coeff_corr)
        target += -0.5 * lambda_spline_corr * square(coeff_corr_eff[m, j] - 2 * coeff_corr_eff[m, j - 1] + coeff_corr_eff[m, j - 2]);
  }
  if (estimate_spline_vcov == 1 && n_coeff_vcov > 2 && lambda_spline_vcov > 0) {
    for (m in 1:M_vcov)
      for (j in 3:n_coeff_vcov)
        target += -0.5 * lambda_spline_vcov * square(coeff_vcov_eff[m, j] - 2 * coeff_vcov_eff[m, j - 1] + coeff_vcov_eff[m, j - 2]);
  }
  if (estimate_spline_cv_mean == 1 && n_coeff_cv_mean > 2 && lambda_spline_cv_mean > 0) {
    for (j in 3:n_coeff_cv_mean)
      target += -0.5 * lambda_spline_cv_mean * square(coeff_cv_mean_eff[j] - 2 * coeff_cv_mean_eff[j - 1] + coeff_cv_mean_eff[j - 2]);
  }
  if (estimate_spline_cv_marker == 1 && n_coeff_cv_marker > 2 && lambda_spline_cv_marker > 0) {
    for (j in 3:n_coeff_cv_marker)
      target += -0.5 * lambda_spline_cv_marker * square(coeff_cv_marker_eff[j] - 2 * coeff_cv_marker_eff[j - 1] + coeff_cv_marker_eff[j - 2]);
  }
  if (estimate_spline_cs_mean == 1 && n_coeff_cs_mean > 2 && lambda_spline_cs_mean > 0) {
    for (j in 3:n_coeff_cs_mean)
      target += -0.5 * lambda_spline_cs_mean * square(coeff_cs_mean_eff[j] - 2 * coeff_cs_mean_eff[j - 1] + coeff_cs_mean_eff[j - 2]);
  }
  if (estimate_spline_cs_marker == 1 && n_coeff_cs_marker > 2 && lambda_spline_cs_marker > 0) {
    for (j in 3:n_coeff_cs_marker)
      target += -0.5 * lambda_spline_cs_marker * square(coeff_cs_marker_eff[j] - 2 * coeff_cs_marker_eff[j - 1] + coeff_cs_marker_eff[j - 2]);
  }


