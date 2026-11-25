#' @section Sample size checks:
#' The model is checked against the number of events it is fitted to. Any
#' problems are reported together in one warning of class
#' `bootlace_sample_size_warning`. The diagnostics are also stored in the
#' returned object as `sample_size` and summarised by the print method.
#'
#' * **Events per variable (EPV).** Events are the number of events for a Cox
#'   model, or the number in the less frequent outcome class for a binary
#'   model. Variables are the estimated coefficients other than the intercept,
#'   so each factor contrast counts. A warning is given when EPV is below 10.
#'   Below that, Peduzzi et al. (1996) found biased regression coefficients
#'   and confidence intervals without proper coverage. [stepwise_reduce()]
#'   uses stricter limits for variable selection.
#' * **Sparse factor levels.** A level of a categorical predictor with no
#'   events (or, for a binary outcome, no non-events). This usually means
#'   separation or an infinite estimate.
#' * **Too few events overall.** Fewer than 20 events in the less frequent
#'   outcome class. This floor is this package's own pragmatic choice, not a
#'   published standard.
#'
#' These are rough rules of thumb. The strength of the effects and the
#' correlation between predictors also affect how accurate the estimates are
#' (Heinze et al., 2018, Section 2.1.3). To silence the warning once you have
#' considered it, wrap the call in
#' `withCallingHandlers(<call>, bootlace_sample_size_warning = function(w) invokeRestart("muffleWarning"))`.
#'
#' @references
#' Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A
#' simulation study of the number of events per variable in logistic
#' regression analysis. *Journal of Clinical Epidemiology*, 49(12),
#' 1373--1379. \doi{10.1016/S0895-4356(96)00236-3}
#'
#' Heinze G, Wallisch C, Dunkler D (2018). Variable selection -- A review and
#' recommendations for the practicing statistician. *Biometrical Journal*,
#' 60(3), 431--449. \doi{10.1002/bimj.201700067}
