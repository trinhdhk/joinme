/**
 * @file functional/parameters/fit.stan
 * @brief Parameters governing fitted association transformations.
 *
 * @details
 * These parameters describe monotone spline increments, piecewise-linear simplexes and optional affine intercepts and slopes inside functional transformations.  They change the shape or affine placement of an association feature, not the longitudinal trajectory itself.
 */
  /**
   * @brief Latent increments for Stan-estimated monotone I-splines.
   *
   * @details These parameters preserve the established I-spline
   * parameterisation and its regularisation-family prior.  Piecewise-linear modes
   * omit them because those modes use the explicit simplexes declared below.
   */
  vector[(tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) ? 0 : n_free_spline_cv] z_spline_cv; // Role: standardised latent value spline current value.
  vector[(tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) ? 0 : n_free_spline_cs] z_spline_cs; // Role: standardised latent value spline current slope.
  matrix[M_corr, (tf_mode_corr == 3 || tf_mode_corr == 7) ? 0 : n_free_spline_corr] z_spline_corr; // Role: standardised latent value spline correlation.
  matrix[M_vcov, (tf_mode_vcov == 3 || tf_mode_vcov == 7) ? 0 : n_free_spline_vcov] z_spline_vcov; // Role: standardised latent value spline covariance.
  vector[(tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) ? 0 : n_free_spline_cv_mean] z_spline_cv_mean; // Role: standardised latent value spline current value mean.
  vector[(tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) ? 0 : n_free_spline_cv_marker] z_spline_cv_marker; // Role: standardised latent value spline current value marker.
  vector[(tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) ? 0 : n_free_spline_cs_mean] z_spline_cs_mean; // Role: standardised latent value spline current slope mean.
  vector[(tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) ? 0 : n_free_spline_cs_marker] z_spline_cs_marker; // Role: standardised latent value spline current slope marker.

  /**
   * @brief Simplex increments for ordered piecewise-linear associations.
   *
   * @details For K ordinates, a K-1 simplex allocates the unit association span
   * between adjacent knots. Cumulative sums in transformed parameters recover
   * the ordered ordinates. An inactive channel has a one-element simplex,
   * which has no free parameter and permits one compiled model to cover every
   * transform choice without changing existing I-spline parameters.
   */
  simplex[(tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) ? n_free_spline_cv : 1] pwlin_simplex_cv; // Role: piecewise-linear simplex current value.
  simplex[(tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) ? n_free_spline_cs : 1] pwlin_simplex_cs; // Role: piecewise-linear simplex current slope.
  array[M_corr] simplex[(tf_mode_corr == 3 || tf_mode_corr == 7) ? n_free_spline_corr : 1] pwlin_simplex_corr; // Role: piecewise-linear simplex correlation.
  array[M_vcov] simplex[(tf_mode_vcov == 3 || tf_mode_vcov == 7) ? n_free_spline_vcov : 1] pwlin_simplex_vcov; // Role: piecewise-linear simplex covariance.
  simplex[(tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) ? n_free_spline_cv_mean : 1] pwlin_simplex_cv_mean; // Role: piecewise-linear simplex current value mean.
  simplex[(tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) ? n_free_spline_cv_marker : 1] pwlin_simplex_cv_marker; // Role: piecewise-linear simplex current value marker.
  simplex[(tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) ? n_free_spline_cs_mean : 1] pwlin_simplex_cs_mean; // Role: piecewise-linear simplex current slope mean.
  simplex[(tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) ? n_free_spline_cs_marker : 1] pwlin_simplex_cs_marker; // Role: piecewise-linear simplex current slope marker.

  /* fit-only affine-shift parameters for functional transforms */
  vector[estimate_iota_intercept_cv] z_iota_intercept_cv; // Role: standardised latent value affine transformation intercept current value.
  vector[estimate_iota_slope_cv] z_iota_slope_cv; // Role: standardised latent value affine transformation slope current value.
  vector[estimate_iota_intercept_cs] z_iota_intercept_cs; // Role: standardised latent value affine transformation intercept current slope.
  vector[estimate_iota_slope_cs] z_iota_slope_cs; // Role: standardised latent value affine transformation slope current slope.
  vector[M_corr * estimate_iota_intercept_corr] z_iota_intercept_corr; // Role: standardised latent value affine transformation intercept correlation.
  vector[M_corr * estimate_iota_slope_corr] z_iota_slope_corr; // Role: standardised latent value affine transformation slope correlation.
  vector[M_vcov * estimate_iota_intercept_vcov] z_iota_intercept_vcov; // Role: standardised latent value affine transformation intercept covariance.
  vector[M_vcov * estimate_iota_slope_vcov] z_iota_slope_vcov; // Role: standardised latent value affine transformation slope covariance.
  vector[estimate_iota_intercept_cv_mean] z_iota_intercept_cv_mean; // Role: standardised latent value affine transformation intercept current value mean.
  vector[estimate_iota_slope_cv_mean] z_iota_slope_cv_mean; // Role: standardised latent value affine transformation slope current value mean.
  vector[estimate_iota_intercept_cv_marker] z_iota_intercept_cv_marker; // Role: standardised latent value affine transformation intercept current value marker.
  vector[estimate_iota_slope_cv_marker] z_iota_slope_cv_marker; // Role: standardised latent value affine transformation slope current value marker.
  vector[estimate_iota_intercept_cs_mean] z_iota_intercept_cs_mean; // Role: standardised latent value affine transformation intercept current slope mean.
  vector[estimate_iota_slope_cs_mean] z_iota_slope_cs_mean; // Role: standardised latent value affine transformation slope current slope mean.
  vector[estimate_iota_intercept_cs_marker] z_iota_intercept_cs_marker; // Role: standardised latent value affine transformation intercept current slope marker.
  vector[estimate_iota_slope_cs_marker] z_iota_slope_cs_marker; // Role: standardised latent value affine transformation slope current slope marker.

