# bootlace 0.1.0

* Fixes from the follow-up scientific review (no issue):
  * `plot_cox_survival(fun = "cumhaz")` returned inverted confidence limits; bounds are now computed on the requested scale.
  * `fp_data()` / `fp_plot()` work with non-treatment contrasts (e.g. `contr.sum`): terms are mapped to their fitted coefficients, and level rows are level-vs-reference comparisons built from the fitted contrasts.
  * `stepwise_reduce()` no longer discards refit warnings: non-convergence makes a comparison untested (listed in `untested`), and separation or possibly infinite coefficients are listed in the new `test_warnings`.
  * Count and seed arguments must be whole numbers within R's integer range.
  * `fit_logistic_emax()` warns with class `bootlace_emax_stopped_warning` when the optimiser stops short of the optimum, separately from a non-identified Hessian; `bootstrap_emax()` keeps both kinds of replicate (so EC50 ridge behaviour stays in the intervals) but counts each reason in `flag_counts`.

* Emax intervals: a simulation study found that bootstrap percentile intervals under-covered the true Emax (86% for a nominal 95% in a well-identified design, 80-83% in a weakly identified one). With the changes below, BCa intervals cover 89-95% and profile intervals 91-99% in both designs:
  * `summary.boot_emax()` and `plot_emax_curve()` now default to BCa (bias-corrected and accelerated) intervals, using jackknife fits that `bootstrap_emax()` stores as `jack_theta`; `ci = "percentile"` gives the old intervals.
  * `bootstrap_emax()` refits each replicate from the automatic multi-start instead of the original estimate, which had understated the spread of EC50 by about 20% in weakly identified designs.
  * New `confint()` method for `fit_logistic_emax()` fits: profile-likelihood intervals by default (with infinite bounds, and a warning, when the data don't determine a parameter in one direction), or Wald intervals.
  * An opt-in coverage study (`BOOTLACE_COVERAGE=true`) checks interval coverage for the logistic, Cox (including patient-level resampling) and Emax bootstraps (no issue).

* Argument validation: bootstrap, simulation and grid counts must be whole numbers; confidence levels, `alpha` and percentiles must lie in (0, 1); `vpc_cox()`'s `pi_probs` must be increasing and `t_max` positive; seeds must be single numbers; `times` and ROC `thresholds` must be finite and in range; `predict.logistic_emax()` rejects negative or infinite exposures (no issue).

* `fp_data()` and `fp_plot()` label interaction rows as ratios of ORs or HRs, and give header rows a joint Wald test over all of the term's coefficients (shown with its df when above 1) instead of the first contrast's p-value. A new `p_df` column records the test's df (no issue).

* `plot_cox_survival()` bands now use each bootstrap replicate's own baseline cumulative hazard (stored by `bootstrap_cox()` as `boot_basehaz`). Holding the baseline fixed made the bands far too wide (for `lung`, age 60, one year: 0.02 to 0.73, versus 0.27 to 0.44 now and 0.28 to 0.45 from `survfit()`). Survival is now 1 before the first event time instead of taking the first step's value (no issue).

* `bootstrap_logistic()`, `bootstrap_cox()` and `bootstrap_emax()` no longer drop every replicate that warns, which conditioned the bootstrap on well-behaved samples. Benign warnings (non-integer weights) are ignored; separation, infinite-coefficient and Emax identifiability warnings keep the replicate and are counted in `n_flagged`; errors, non-convergence, other warnings and `NA` coefficients drop it, with reasons in `failures`. Results also record `n_requested` (no issue).

* bootlace now requires ggplot2 (>= 4.0.0). `roc_logistic_emax()` and `plot_cox_forest()` no longer use arguments and geoms that ggplot2 4.0 deprecates (no issue).
* `plot_logistic_curves()` with `colour_by`: categorical predictors with `obs_by_group = TRUE` keep their per-level marker colours (the observed squares' fill scale used to replace them), and continuous predictors no longer get an unused fill scale that ggplot2 4.0 warns about (no issue).

* `bootstrap_logistic()` and `bootstrap_cox()` refit from the model's own design, so the analysis rows, case weights, offsets and strata can't change between replicates. `bootstrap_cox()` resamples whole subjects when the fit has an `id` or cluster, and rejects counting-process data without one (no issue).
* `stepwise_reduce()` refits candidate models on the full model's rows with its weights and offset, and computes the GLM likelihood-ratio statistic as the deviance difference, as `drop1()` does. This gives the same result as before except with non-integer case weights, where `logLik()` of a binomial GLM is not a true log-likelihood (no issue).
* `roc_logistic()` and `roc_logistic_emax()` report the exact empirical AUC (ties count one half) and take the Youden-optimal threshold from every observed score; `thresholds` now only controls what is drawn. `roc_logistic()` checks that `outcome` lines up with each model's fitted values (no issue).
* `fit_logistic_emax()` validates its inputs, refuses designs with fewer than 3 distinct exposures, tries several EC50 starting values, and reports `converged = TRUE` and standard errors only when the Hessian is positive definite and the optimum is reached. It gives a detailed warning (class `bootlace_emax_ec50_warning`) when EC50 lies above the highest exposure or below a tenth of the lowest non-zero exposure. `bootstrap_emax()` keeps replicates with such EC50 values instead of dropping them, and summarises them in one detailed warning, including whether they set the percentile interval bounds (no issue).
* `fp_data()` and `fp_plot()` accept only Cox and logit-link binomial models, the ones whose exponentiated coefficients are hazard or odds ratios (no issue).
* `vpc_cox()` checks that `data` holds the fitted records in the fitted order, and bootstraps whole clusters when the fit has one (no issue).

* `bootstrap_logistic()`, `bootstrap_cox()`, `fit_logistic_emax()`, `stepwise_reduce()` and `vpc_cox()` check the number of events against the model. They warn (class `bootlace_sample_size_warning`) when events per variable is below 10, at most 25 for variable selection, when a factor level has no events or non-events, when there are fewer than 20 events, or (Emax) when there are fewer than 4 distinct exposure levels. The diagnostics are stored as `sample_size` and shown by the print methods (no issue).

* Initial package development for bootstrap-based inference and response-curve
  visualization for binomial logistic regression models.
* `compare_logistic_roc()` adds overlaid ROC-curve comparison plots for multiple fitted binomial GLMs (no issue).
* `vpc_cox()` adds a simulation-based visual predictive check for `coxph()` fits, comparing observed Kaplan-Meier curves with prediction intervals of simulated ones, optionally stratified by a grouping column and with bootstrap parameter uncertainty (no issue).
