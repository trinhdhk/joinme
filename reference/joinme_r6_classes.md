# JoiNMe R6 Classes

R6 classes for JoiNMe fit and prediction results with mutable state.

**JoiNMeFit**: holds fitted model, Stan data, and formulas. Public
fields: fit, stan_data, formulaLong, formulaEvent, formulaVCov, config,
call, tmax, dataLong, dataEvent. Methods: initialize, cache_get,
cache_set.

**JoiNMeDynPred**: holds dynamic predictions and metadata. Public
fields: predictions, quantiles, draws, data, metadata, call, tmax,
n_samples. Methods: initialize, cache_get, cache_set.

**SummaryJoiNMeFit**: holds cached summary tables and diagnostics.
Public fields: tables, diagnostics, metadata. Methods: initialize.

**SummaryJoiNMeDynPred**: holds summary tables for predictions. Public
fields: tables, metadata. Methods: initialize.
