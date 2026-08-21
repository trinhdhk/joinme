/**
 * @file functional/transformed_parameters/fit.stan
 * @brief Construct fitted transformation shapes and affine shifts.
 *
 * @details
 * Simplex or softmax increments produce unit-span monotone curves.  Fitted functional intercepts and slopes are then restored from their component-specific prior positions.  Unit-span shapes keep transformation scale separate from the hazard association coefficient.
 */
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
    vector[n_iota_prior] iota_prior_effect; // scientific affine shifts extracted from the common transformed prior programme
    int iota_position = 1; // next position within the functional component's local slice

    if (n_iota_prior > 0)
      iota_prior_effect = regression_prior_effect[prior_start_iota : prior_start_iota + n_iota_prior - 1];
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


