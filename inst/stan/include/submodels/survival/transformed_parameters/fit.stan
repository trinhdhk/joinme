/**
 * @file survival/transformed_parameters/fit.stan
 * @brief Document the survival transformed-parameter stage.
 *
 * @details
 * The survival submodel needs no additional derived parameter after the common coefficient programme has restored event-covariate slopes.  Baseline-hazard coefficients remain on their sampled scale.  Event-time calculations occur in the joint likelihood because they require longitudinal association features.
 */



