#' bootlace: Bootstrap Inference and Plots for Logistic, Emax and Cox Models
#'
#' Provides case-resampling bootstrap for fitted binomial GLMs, logistic Emax
#' models and Cox proportional hazards models, with companion plots and model
#' diagnostics.
#'
#' The main entry points are:
#'
#' **Logistic regression**
#' - [bootstrap_logistic()]: bootstrap a fitted binomial `glm`.
#' - [plot_logistic_curves()]: response-curve plots from the result.
#' - [roc_logistic()], [compare_logistic_roc()]: ROC curves for one or more
#'   models.
#'
#' **Logistic Emax models**
#' - [fit_logistic_emax()], [bootstrap_emax()]: fit and bootstrap the model.
#' - [plot_emax_curve()]: response-curve plot from the result.
#' - [roc_logistic_emax()]: ROC curve for the fit.
#'
#' **Cox proportional hazards models**
#' - [bootstrap_cox()]: bootstrap a fitted `coxph` model.
#' - [plot_cox_forest()], [plot_cox_survival()]: hazard-ratio forest plot and
#'   adjusted survival curves from the result.
#' - [vpc_cox()]: visual predictive check against observed Kaplan-Meier
#'   curves.
#'
#' **Binomial GLMs and Cox models**
#' - [fp_data()], [fp_plot()]: build and render a forest plot.
#' - [stepwise_reduce()]: backwards stepwise reduction by likelihood ratio
#'   test.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats AIC BIC coef formula glm model.frame plogis pnorm qnorm
#'   quantile rbinom rnorm setNames
## usethis namespace: end
NULL

## Suppress R CMD CHECK NOTEs for ggplot2 aesthetic names used inside aes()
utils::globalVariables(
  c(
    "x", "y", "prob", "ci_lo", "ci_hi", "obs", "grp", "D", "M", "model",
    "fpr", "sensitivity", "threshold", "specificity",
    "term", "hr", "time", "surv", "curve",
    "lo", "hi", "median", "group",
    # refit columns and variables evaluated by glm()/coxph() via NSE
    ".bootlace_w", ".bootlace_off", "W"
  )
)
