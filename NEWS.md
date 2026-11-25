# bootlace 0.1.0

* Initial package development for bootstrap-based inference and response-curve
  visualization for binomial logistic regression models.
* `compare_logistic_roc()` adds overlaid ROC-curve comparison plots for multiple fitted binomial GLMs (no issue).
* `vpc_cox()` adds a simulation-based visual predictive check for `coxph()` fits, comparing observed Kaplan-Meier curves with prediction intervals of simulated ones, optionally stratified by a grouping column and with bootstrap parameter uncertainty (no issue).
