/**
 * @file include/submodels/latent_class/data/fit.stan
 * @brief Data declarations used only by latent-progress mixture fitting.
 *
 * @details The master mixture programme first includes the five scientific
 * data fragments and `regression_prior_data.stan`. Keeping these class fields
 * separate means that ordinary fitting neither parses mixture dimensions nor
 * receives mixture arrays.
 */

/**
 * @brief Marginalised latent-progress mixture controls.
 *
 * @details This declaration module is included only by `joinme_mix()`.
 * The ordinary fitting programme does not parse or allocate these fields.
 * The mixture supplies one coordinate segment per selected random-effect
 * level.  Starts are one-based when active and zero when inactive.
 */
int<lower=0, upper=1> use_mixture; // whether latent-progress classes are fitted
int<lower=1> n_classes; // number of shared latent classes
int<lower=0> K_mix; // total selected latent-coordinate count
int<lower=0, upper=2> shrinkage; // latent-progress component family: 0 Student-t(6), 1 Laplace, 2 Normal
int<lower=0, upper=1> mix_subject; // whether subject effects receive a latent-class distribution
int<lower=0, upper=1> mix_marker; // whether marker effects receive a latent-class distribution
int<lower=0, upper=1> mix_covariance; // whether covariance effects receive a latent-class distribution

int<lower=0> mix_dim_subject; // selected subject-effect coordinates
array[mix_dim_subject] int<lower=1, upper=R_id> mix_idx_subject; // source indices of selected subject effects
int<lower=0, upper=K_mix> mix_start_subject; // first packed subject coordinate

int<lower=0> mix_dim_marker; // selected marker-effect coordinates
array[mix_dim_marker] int<lower=1, upper=R_mk> mix_idx_marker; // source indices of selected marker effects
int<lower=0, upper=K_mix> mix_start_marker; // first packed marker coordinate

int<lower=0> mix_dim_covariance; // selected covariance-regression coordinates
array[mix_dim_covariance] int<
  lower=1,
  upper=(indep_idmarker_cov == 1
    ? Q_idm
    : (Q_idm * (Q_idm + 1)) %/% 2)
> mix_idx_covariance; // source indices of selected covariance effects
int<lower=0, upper=K_mix> mix_start_covariance; // first packed covariance coordinate

int<lower=0, upper=2> mix_ordering; // label rule: none, random intercept, or baseline probability
int<lower=0, upper=K_mix> mix_ordered_location_coordinate; // sole packed intercept coordinate ordered across classes
vector<lower=0>[n_classes] mix_probability_prior; // Dirichlet concentration from jm_prior()
int<lower=0> P_class_subject; // concatenated subject class predictors
matrix[n_id, P_class_subject] X_class_subject; // subject class-regression design
array[n_classes] int<lower=1, upper=max(1, P_class_subject)> class_term_start_subject; // first subject coefficient used by each class
array[n_classes] int<lower=0, upper=P_class_subject> class_term_count_subject; // subject coefficient count used by each class
int<lower=0> P_class_marker; // concatenated marker class predictors
matrix[D, P_class_marker] X_class_marker; // marker class-regression design
array[n_classes] int<lower=1, upper=max(1, P_class_marker)> class_term_start_marker; // first marker coefficient used by each class
array[n_classes] int<lower=0, upper=P_class_marker> class_term_count_marker; // marker coefficient count used by each class

/**
 * Block-specific prior for the concatenated class-membership regressions.
 * Subject-domain coefficients precede marker-domain coefficients. Family
 * codes match the common prior module: 1 Student-t, 2 Normal, 3 Laplace and
 * 4 regularised horseshoe.
 */
int<lower=1, upper=4> prior_class_regression_family; // selected family for formulaClass coefficients
vector[P_class_subject + P_class_marker] prior_class_regression_mu; // coefficient-specific locations in subject-then-marker order
vector<lower=0>[P_class_subject + P_class_marker] prior_class_regression_scale; // coefficient-specific scales in the same order
real<lower=0> prior_class_regression_df; // fixed Student-t or horseshoe-local degrees of freedom
real<lower=0> prior_class_regression_global_df; // horseshoe global-scale degrees of freedom
real<lower=0> prior_class_regression_global_scale; // horseshoe global scale
real<lower=0> prior_class_regression_slab_df; // horseshoe finite-slab degrees of freedom
real<lower=0> prior_class_regression_slab_scale; // horseshoe finite-slab scale
