# Parameter-stage helpers

Ordinary parameters are declared by scientific submodels.
`regression_prior_parameters.stan` contains only shrinkage auxiliaries shared
by several coefficient families, while
`dynamic_prediction_subject_effects.stan` contains new-subject coordinates
coupled across prediction targets. Mixture-only parameters remain separate so
they do not inflate the ordinary Stan programme.
