/**
 * @file include/etc/parameters/dynamic_prediction_subject_effects.stan
 * @brief Sample new-subject latent effects conditionally for every fitted draw.
 *
 * @details
 * The single array retains the established draw-major storage expected by the
 * threaded conditional kernel.  Its internal offsets distinguish subject,
 * marker and covariance latent coordinates.
 */
  /* Subject-specific random effects (standard normal to be scaled) */
  array[n_draws] vector[n_random_id] z_u; // id-level latent normals per draw
  array[n_draws, n_marker_types] vector[n_random_marker] z_v; // marker latents
  array[n_draws, n_marker_types] vector[n_random_marker_id] z_w_lat; // marker-id latents
  array[n_draws] vector[num_unique_cov_entries] z_L; // latent normals for covariance regression

