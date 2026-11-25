# bootlace 0.1.0

* `bootstrap_logistic()`, `bootstrap_cox()`, `fit_logistic_emax()`, `stepwise_reduce()` and `vpc_cox()` check the number of events against the model. They warn (class `bootlace_sample_size_warning`) when events per variable is below 10, at most 25 for variable selection, when a factor level has no events or non-events, when there are fewer than 20 events, or (Emax) when there are fewer than 3 non-zero exposure levels. The diagnostics are stored as `sample_size` and shown by the print methods (no issue).

* Initial package development for bootstrap-based inference and response-curve
  visualization for binomial logistic regression models.
* `compare_logistic_roc()` adds overlaid ROC-curve comparison plots for multiple fitted binomial GLMs (no issue).
* `vpc_cox()` adds a simulation-based visual predictive check for `coxph()` fits, comparing observed Kaplan-Meier curves with prediction intervals of simulated ones, optionally stratified by a grouping column and with bootstrap parameter uncertainty (no issue).
