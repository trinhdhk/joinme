/**
 * @file assoc/data/dynamic_prediction.stan
 * @brief Declare posterior association coefficients and active channels.
 *
 * @details
 * Each retained draw supplies scalar, correlation and covariance association coefficients.  Flags retain the fitted association specification when future hazards are evaluated.
 */
  array[n_draws] real coeff_assoc_cv_total;     // association coeffs: total CV
  array[n_draws] real coeff_assoc_cs_total;     // association coeffs: total CS
  array[n_draws] real coeff_assoc_cv_mean;      // association coeffs: mean CV
  array[n_draws] real coeff_assoc_cs_mean;      // association coeffs: mean CS
  array[n_draws] real coeff_assoc_cv_marker;    // association coeffs: marker CV
  array[n_draws] real coeff_assoc_cs_marker;    // association coeffs: marker CS
  array[n_draws] vector[(n_random_marker_id * (n_random_marker_id - 1)) %/% 2] coeff_assoc_corr; // corr coeffs (off-diagonal correlations)
  array[n_draws] vector[num_unique_cov_entries] coeff_assoc_vcov; // vcov coeffs for lower-triangular L entries

  int<lower=0, upper=1> flag_assoc_cv_total;    // include total CV association
  int<lower=0, upper=1> flag_assoc_cv_mean;     // include mean CV association
  int<lower=0, upper=1> flag_assoc_cv_marker;   // include marker CV association
  int<lower=0, upper=1> flag_assoc_cs_total;    // include total CS association
  int<lower=0, upper=1> flag_assoc_cs_mean;     // include mean CS association
  int<lower=0, upper=1> flag_assoc_cs_marker;   // include marker CS association
  int<lower=0, upper=1> flag_assoc_corr;        // include corr association
  int<lower=0, upper=1> flag_assoc_vcov;        // include vcov association
  int<lower=0> M_corr_tf;                       // number of corr transform components
  int<lower=0> M_vcov_tf;                       // number of vcov transform components



