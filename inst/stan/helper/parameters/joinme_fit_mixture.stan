/**
 * @file joinme_fit_mixture.stan
 * @brief Parameters belonging exclusively to the latent-progress mixture.
 *
 * @details Component allocations are marginalised, so no discrete parameter
 * is declared. The simplex, locations and positive scales are shared by all
 * selected clustering blocks through the common coordinate layout.
 */

/**
 * @brief Parameters of the marginalised latent-progress mixture.
 *
 * @details `simplex[1]` has no free coordinate and a matrix with zero columns
 * has no elements, so the ordinary non-mixture model retains its original
 * posterior dimension.
 */
array[mix_ordering == 2 ? 1 : 0] positive_ordered[n_clusters - 1] mix_probability_ordered_gap; // positive gaps giving ordered log-odds below the final class
array[mix_ordering == 2 ? 0 : 1] simplex[n_clusters] mix_probability_unordered; // unrestricted baseline class probabilities
array[mix_ordered_location_coordinate > 0 ? 1 : 0] ordered[n_clusters] mix_location_ordered; // ordered class locations for the selected random intercept only
matrix[n_clusters, K_mix - (mix_ordered_location_coordinate > 0 ? 1 : 0)] mix_location_unordered; // unrestricted slope and other selected locations
matrix<lower=1e-8>[n_clusters, K_mix] mix_scale; // positive component scales by class and coordinate
vector[P_class_subject] mix_class_coefficient_subject; // concatenated subject-domain class coefficients
vector[P_class_marker] mix_class_coefficient_marker; // concatenated marker-domain class coefficients
