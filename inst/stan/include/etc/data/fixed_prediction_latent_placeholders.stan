/**
 * @file include/etc/data/fixed_prediction_latent_placeholders.stan
 * @brief Declare neutral latent arrays for the fixed fitted-effect programme.
 *
 * @details
 * The common prediction calculation names the new-subject latent arrays even
 * though the fitted-effect route selects realised posterior effects instead.
 * Declaring zero-valued arrays as data permits that common scientific
 * calculation to be reused without introducing parameters or a model block.
 */
array[n_draws] vector[n_random_id] z_u; // unused neutral subject latents
array[n_draws, n_marker_types] vector[n_random_marker] z_v; // unused neutral marker latents
array[n_draws, n_marker_types] vector[n_random_marker_id] z_w_lat; // unused neutral marker-by-subject latents
array[n_draws] vector[num_unique_cov_entries] z_L; // unused neutral covariance-regression latents
