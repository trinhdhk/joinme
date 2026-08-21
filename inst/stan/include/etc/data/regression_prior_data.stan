/**
 * @file include/etc/data/regression_prior_data.stan
 * @brief Cross-submodel data for the common coefficient-prior programme.
 *
 * @details
 * R has already classified each fitted coefficient by scientific submodel and
 * by intercept or slope role.  These arrays retain that classification after
 * packing, including coefficient-wise families and regularised-horseshoe group
 * indices. This helper owns the information because a single transformation
 * is shared by longitudinal, survival, association and functional
 * coefficients, and therefore does not belong to one scientific submodel.
 */
  /**
   * Block-specific coefficient prior declarations.
   * Family codes are 1 Student-t, 2 Normal, 3 Laplace and 4 regularised
   * horseshoe. Marker and marker-weight declarations contain no ordinary
   * location or scale data: their transformations construct zero and one
   * explicitly, so only the family and its tail/shrinkage hyperparameters can
   * vary.
   */
  /* Common regression-prior programme: component roles are resolved in R. */
  int<lower=0> n_regression_prior; // total number of longitudinal, event, covariance, association, functional and distributional coefficients
  array[n_regression_prior] int<lower=1, upper=4> prior_regression_family; // coefficient-wise family: Student-t, Normal, Laplace or regularised horseshoe
  vector[n_regression_prior] prior_regression_mu; // coefficient-wise scientific prior locations
  vector<lower=0>[n_regression_prior] prior_regression_scale; // coefficient-wise ordinary prior scales
  vector<lower=0>[n_regression_prior] prior_regression_df; // Student-t degrees of freedom or horseshoe local degrees of freedom
  array[n_regression_prior] int<lower=0> prior_regression_horseshoe_local_index; // local-scale index, zero for coefficients outside a horseshoe group
  array[n_regression_prior] int<lower=0> prior_regression_horseshoe_group_index; // global-scale group index, zero outside a horseshoe group
  int<lower=0> n_regression_horseshoe_local; // number of coefficient-specific horseshoe local scales
  vector<lower=0>[n_regression_horseshoe_local] prior_regression_horseshoe_local_df; // half-Student-t degrees of freedom for local scales
  int<lower=0> n_regression_horseshoe_group; // number of component-role groups using a horseshoe
  vector<lower=0>[n_regression_horseshoe_group] prior_regression_horseshoe_global_df; // half-Student-t degrees of freedom for group-global scales
  vector<lower=0>[n_regression_horseshoe_group] prior_regression_horseshoe_global_scale; // fixed global shrinkage scales
  vector<lower=0>[n_regression_horseshoe_group] prior_regression_horseshoe_slab_df; // finite-slab degrees of freedom by group
  vector<lower=0>[n_regression_horseshoe_group] prior_regression_horseshoe_slab_scale; // finite-slab scales by group
  int<lower=0> prior_start_beta; // first longitudinal population coefficient in the packed prior vector
  int<lower=0> prior_start_alpha; // first association coefficient in the packed prior vector
  int<lower=0> prior_start_iota; // first functional affine-shift coefficient in the packed prior vector
  int<lower=0> prior_start_vcov_sd; // first standard-deviation covariance-regression coefficient
  int<lower=0> prior_start_vcov_corr; // first partial-correlation covariance-regression coefficient
  int<lower=0> prior_start_survival; // first event-covariate coefficient
  int<lower=0> prior_start_sigma; // first sigma-regression coefficient
  int<lower=0> prior_start_nu; // first degrees-of-freedom regression coefficient
  int<lower=0> prior_start_phi; // first dispersion-regression coefficient
  int<lower=0> prior_start_distributional_alpha; // first skewness-regression coefficient
  int<lower=0> prior_start_kappa; // first beta-sample-size regression coefficient
  int<lower=0> prior_start_tau; // first quantile-regression coefficient


