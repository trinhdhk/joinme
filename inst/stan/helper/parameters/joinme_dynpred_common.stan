  /* Subject-specific random effects (standard normal to be scaled) */
  array[n_draws] vector[n_random_id] z_u; // id-level latent normals per draw
  array[n_draws, n_marker_types] vector[n_random_marker] z_v; // marker latents
  array[n_draws, n_marker_types] vector[n_random_marker_id] z_w_lat; // marker-id latents
  vector[n_draws] z_L; // latent normals for covariance regression
