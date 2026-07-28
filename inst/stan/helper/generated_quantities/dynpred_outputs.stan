/**
 * @file dynpred_outputs.stan
 * @brief Compatibility include for ordinary dynamic-prediction outputs.
 *
 * @details New entry points include the declaration and calculation modules
 * directly. This adapter preserves the established helper path for downstream
 * code which includes the ordinary output bundle.
 */
#include dynpred_output_declarations.stan
#include dynpred_output_calculations.stan
