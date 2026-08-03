/**
 * @file joinme_fit_mixture.stan
 * @brief Parameters belonging exclusively to the latent-progress mixture.
 *
 * @details Component allocations are marginalised, so no discrete parameter
 * is declared. The simplex, locations and positive scales are shared by all
 * selected class-specific blocks through the common coordinate layout.
 *
 * Locations and scales live on the standardised latent coordinates `z_u`,
 * `z_v`, and `z_L`. They are not standard deviations on the response scale.
 * The established transformations subsequently map `z_u` through `L_u`,
 * `z_v` through `L_v`, and `z_L` through the non-negative loading `lambda_L`.
 * `simulate_joinme_mix()` follows that same order: it first draws a
 * component-conditional standardised coordinate and only then applies the
 * ordinary random-effect transformation.
 *
 * Consequently, `mix_scale` and the ordinary scales `tau_u`, `tau_v`, or
 * `lambda_L` can trade off: the likelihood often learns their product much
 * more directly than either factor. This is statistical weak scale
 * identification, not a difference between simulation and fitting. It may
 * produce a long posterior ridge, low E-BFMI and high autocorrelation when the
 * data contain little information that separates the two factors. The direct
 * Exponential(1) prior on `mix_scale`, requested by the model specification,
 * is stated in `fit_mixture_priors.stan`; it regularises but does not remove
 * that likelihood ridge. `diagnosis()` reports the resulting chain-specific
 * energy behaviour and posterior scale correlations explicitly.
 */

/**
 * @brief Parameters of the marginalised latent-progress mixture.
 *
 * @details `simplex[1]` has no free coordinate and a matrix with zero columns
 * has no elements, so the ordinary non-mixture model retains its original
 * posterior dimension.
 */
array[mix_ordering == 2 ? 1 : 0] positive_ordered[n_classes - 1] mix_probability_ordered_gap; // positive gaps giving ordered log-odds below the final class
array[mix_ordering == 2 ? 0 : 1] simplex[n_classes] mix_probability_unordered; // unrestricted baseline class probabilities
array[mix_ordered_location_coordinate > 0 ? 1 : 0] ordered[n_classes] mix_location_ordered; // ordered class locations for the selected random intercept only
matrix[n_classes, K_mix - (mix_ordered_location_coordinate > 0 ? 1 : 0)] mix_location_unordered; // unrestricted slope and other selected locations
matrix<lower=1e-8>[n_classes, K_mix] mix_scale; // positive within-class component scales by class and packed coordinate
vector[P_class_subject] mix_class_coefficient_subject; // concatenated subject-domain class coefficients
vector[P_class_marker] mix_class_coefficient_marker; // concatenated marker-domain class coefficients
