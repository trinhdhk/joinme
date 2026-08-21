# Transformed-data helpers

Deterministic preparation belongs to a scientific submodel whenever possible.
`dynamic_prediction_draw_indices.stan` constructs the common retained-draw
sequence required by the joint threaded prediction calculation. All other
parameter-free preparation remains in its scientific submodel.
