/**
 * @file marker_weight/data/fit.stan
 * @brief Data declarations for hierarchical marker weighting.
 *
 * @details
 * Marker-weight data vectors are known offsets declared inside
 * priors$marker_weights in R. Set indices state whether active
 * weighted association terms share a weight vector.  Estimated weights add one
 * fitted common location per set, together with one centred, unit-scale
 * departure per marker. The common location receives the ordinary regression prior declared
 * by marker_weights$intercept in R. The prior data below therefore describe
 * only the departure family and
 * its tail or regularised-horseshoe hyperparameters. For Student-t
 * departures, estimate_marker_weight_df distinguishes a fixed
 * prior_student_t() declaration from the family name `student_t`, which learns
 * degrees of freedom separately for every active weight set. The public R
 * families `constant` and `none` set estimate_marker_weights to zero; the
 * offset vectors are then already the complete effective weights. No extra
 * Stan family code is needed because every marker-weight parameter dimension
 * is zero in that case.
 *
 * The longitudinal data fragment must be included first because it declares D.
 */
  vector[D] marker_weights_cv_total;   // known offset for total current-value association
  vector[D] marker_weights_cs_total;   // known offset for total current-slope association
  vector[D] marker_weights_cv_marker;  // known offset for marker current-value association
  vector[D] marker_weights_cs_marker;  // known offset for marker current-slope association
  int<lower=0, upper=1> marker_weight_sets_shared; // 1 if weighted association terms share one marker-weight structure
  int<lower=1, upper=4> n_marker_weight_sets;  // number of distinct shared or term-specific marker-weight sets
  int<lower=0, upper=4> marker_weight_set_cv_total;  // set index for cv_total weights, 0 when inactive
  int<lower=0, upper=4> marker_weight_set_cs_total;  // set index for cs_total weights, 0 when inactive
  int<lower=0, upper=4> marker_weight_set_cv_marker; // set index for cv_marker weights, 0 when inactive
  int<lower=0, upper=4> marker_weight_set_cs_marker; // set index for cs_marker weights, 0 when inactive
  int<lower=0, upper=1> estimate_marker_weights; // 1 fits a common location and direct unit-scale signed departures around offsets; 0 uses offsets exactly
  int<lower=0, upper=1> use_marker_weight_assoc; // 1 when marker-weighted assoc terms are active
  int<lower=0, upper=4> n_marker_weight_means; // number of fitted common marker-weight locations; zero for fixed or inactive weights
  int<lower=1, upper=4> prior_marker_weight_family; // family of standardised marker-specific weight departures
  int<lower=0, upper=1> estimate_marker_weight_df; // 1 learns one Student-t degrees of freedom per active weight set; 0 uses an explicitly declared fixed value or a non-Student family
  real<lower=0> prior_marker_weight_df; // fixed Student-t degrees of freedom when declared, or local-scale degrees of freedom for the horseshoe family
  real<lower=0> prior_marker_weight_global_df; // unit-default horseshoe global degrees of freedom for marker weights
  real<lower=0> prior_marker_weight_global_scale; // configured horseshoe global shrinkage scale for marker-weight departures
  real<lower=0> prior_marker_weight_slab_df; // fixed horseshoe slab degrees of freedom for marker weights
  real<lower=0> prior_marker_weight_slab_scale; // configured finite-slab scale for marker-weight departures
