# bootlace 0.1.0

* Initial CRAN submission.

* Bootstrap inference for binomial GLMs (`bootstrap_logistic()`), logistic
  Emax models (`fit_logistic_emax()`, `bootstrap_emax()`) and Cox models
  (`bootstrap_cox()`), with percentile intervals; Emax fits also get BCa and
  profile-likelihood intervals (`summary.boot_emax()`, `confint()`), checked
  by an opt-in coverage simulation study.

* Plots: response curves with bootstrap ribbons (`plot_logistic_curves()`,
  `plot_emax_curve()`), forest plots (`fp_plot()`, `fp_data()`,
  `plot_cox_forest()`), adjusted survival curves (`plot_cox_survival()`), ROC
  curves (`roc_logistic()`, `compare_logistic_roc()`, `roc_logistic_emax()`)
  and visual predictive checks for Cox models (`vpc_cox()`).

* Backwards stepwise reduction by likelihood ratio test for GLMs and Cox
  models (`stepwise_reduce()`).

* Events-per-variable checks warn when a model has too few events for its
  number of parameters, or a stepwise search too few for its candidate set.
