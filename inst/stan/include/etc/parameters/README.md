# Parameter-stage helpers

Ordinary parameters are declared by scientific submodels.
`regression_prior_parameters.stan` contains only regularisation auxiliaries shared
by several coefficient families, while
`dynamic_prediction_subject_effects.stan` contains new-subject coordinates
coupled across prediction targets.
