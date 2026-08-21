# Generated-quantity helpers

Ordinary and latent-class output declarations are component-owned. This
directory retains only joint dynamic-prediction calculations whose draw-level
dependence spans several submodels.

`dynpred_output_calculations.stan` selects either newly conditioned effects or
paired realised fitted effects before evaluating the common longitudinal and
survival quantities. The latter selection is used only by the parameter-free
fitted-effect master programme.
