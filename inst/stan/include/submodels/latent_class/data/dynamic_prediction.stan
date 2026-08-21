/**
 * @file include/submodels/latent_class/data/dynamic_prediction.stan
 * @brief Draw-specific mixture data used only by latent-class prediction.
 *
 * @details Each retained fitted-model draw contributes its component
 * probabilities, locations and scales. The layout indices preserve the
 * allocation domains used during fitting.
 */

/**
 * @brief Fitted latent-progress distribution used for dynamic conditioning.
 *
 * @details
 * Dynamic prediction samples latent effects for the subject whose history
 * is being conditioned upon.  When the parent model was fitted by
 * `joinme_mix()`, selected coordinates must therefore use the fitted
 * component distribution rather than silently returning to a standard
 * Normal prior.  The component probability, location and scale are supplied
 * separately for every retained fitted-model draw.
 *
 * This declaration module is loaded only by latent-class dynamic
 * prediction, so an ordinary prediction programme has no mixture-shaped
 * data and cannot accidentally allocate those arrays.
 */
int<lower=0, upper=1> use_dynamic_mixture; // whether the parent fit used latent progress classes
int<lower=1> dynamic_n_classes; // common number of fitted latent classes
int<lower=1> dynamic_mix_dimension; // selected-coordinate width, with one inert column when inactive
int<lower=0, upper=2> dynamic_mix_family; // component family: Student-t(6), Laplace, or Normal
matrix<lower=0>[n_draws, dynamic_n_classes]
  dynamic_mix_probability_subject; // formula-adjusted class probabilities for the new subject
array[n_draws] matrix<lower=0>[
  n_marker_types,
  dynamic_n_classes
] dynamic_mix_probability_marker; // formula-adjusted class probabilities for each fitted marker
array[n_draws] matrix[
  dynamic_n_classes,
  dynamic_mix_dimension
] dynamic_mix_location; // component locations paired with each retained fit draw
array[n_draws] matrix<lower=1e-8>[
  dynamic_n_classes,
  dynamic_mix_dimension
] dynamic_mix_scale; // positive component scales paired with each retained fit draw

// A subject allocation is common to the selected shared-subject and
// covariance-regression coordinates.
int<lower=0, upper=1> dynamic_mix_subject; // whether shared-subject effects receive a latent-class distribution
int<lower=0> dynamic_mix_dim_subject; // number of selected shared-subject coordinates
array[dynamic_mix_dim_subject] int<
  lower=1,
  upper=n_random_id
> dynamic_mix_idx_subject; // selected positions in the shared-subject effect
int<lower=0, upper=dynamic_mix_dimension>
  dynamic_mix_start_subject; // first shared-subject column in the common mixture layout

int<lower=0, upper=1> dynamic_mix_covariance; // whether covariance-regression effects receive a latent-class distribution
int<lower=0> dynamic_mix_dim_covariance; // number of selected covariance-regression coordinates
array[dynamic_mix_dim_covariance] int<
  lower=1,
  upper=num_unique_cov_entries
> dynamic_mix_idx_covariance; // selected positions in the covariance-regression effect
int<lower=0, upper=dynamic_mix_dimension>
  dynamic_mix_start_covariance; // first covariance-regression column in the common mixture layout

// A marker allocation is common to its selected marker effects.
int<lower=0, upper=1> dynamic_mix_marker; // whether shared-marker effects receive a latent-class distribution
int<lower=0> dynamic_mix_dim_marker; // number of selected shared-marker coordinates
array[dynamic_mix_dim_marker] int<
  lower=1,
  upper=n_random_marker
> dynamic_mix_idx_marker; // selected positions in the shared-marker effect
int<lower=0, upper=dynamic_mix_dimension>
  dynamic_mix_start_marker; // first shared-marker column in the common mixture layout

