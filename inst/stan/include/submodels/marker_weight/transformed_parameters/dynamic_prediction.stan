/**
 * @file marker_weight/transformed_parameters/dynamic_prediction.stan
 * @brief Explain the marker-weight dynamic-prediction transformed parameters stage.
 *
 * @details
 * Effective weights enter the conditional association calculation directly
 * from data. They already contain the marker-weight location and direct
 * unit-scale departure reconstructed during fitting, so no further
 * transformation is declared in this programme.
 */
