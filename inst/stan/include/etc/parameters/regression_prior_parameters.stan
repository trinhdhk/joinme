/**
 * @file include/etc/parameters/regression_prior_parameters.stan
 * @brief Parameters shared by the common coefficient-prior programme.
 *
 * @details
 * Coefficient-local horseshoe scales and component-role global and slab scales are stored once because the prior interpreter spans several scientific submodels without giving them a common shrinkage group.
 */
  /**
   * Regularised-horseshoe auxiliaries for the common regression-prior
   * programme.  There is one local scale for each coefficient assigned to a
   * horseshoe and one global/slab pair for each scientific component-role
   * group.  The integer maps in the data block connect coefficients to these
   * vectors without introducing unused parameters for other prior families.
   */
  vector<lower=0>[n_regression_horseshoe_local] horseshoe_local_regression; // coefficient-specific local shrinkage scales
  vector<lower=0>[n_regression_horseshoe_group] horseshoe_global_regression; // component-role global shrinkage scales
  vector<lower=0>[n_regression_horseshoe_group] horseshoe_slab_regression; // finite-slab variance multipliers by component-role group

