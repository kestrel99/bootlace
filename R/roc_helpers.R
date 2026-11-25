# Shared ROC calculations for roc_logistic() and roc_logistic_emax().
#
# AUC and the Youden-optimal cut-off always come from the exact empirical ROC
# curve (every observed score is a cut-off). User-supplied thresholds only
# choose the points that are drawn.

# Empirical ROC at the given cut-offs (classify positive when score >= cut-off).
# With `thresholds = NULL` every observed score is used, plus Inf for (0, 0).
.roc_curve <- function(y_obs, y_pred, thresholds = NULL) {
  if (is.null(thresholds))
    thresholds <- c(Inf, sort(unique(y_pred), decreasing = TRUE))
  n1 <- sum(y_obs == 1)
  n0 <- sum(y_obs == 0)
  tp <- vapply(thresholds, function(t) sum(y_pred >= t & y_obs == 1), numeric(1))
  fp <- vapply(thresholds, function(t) sum(y_pred >= t & y_obs == 0), numeric(1))
  df <- data.frame(
    threshold   = thresholds,
    sensitivity = tp / n1,
    specificity = 1 - fp / n0
  )
  df$fpr    <- 1 - df$specificity
  df        <- df[order(df$fpr, df$sensitivity), ]
  df$youden <- df$sensitivity + df$specificity - 1
  rownames(df) <- NULL
  df
}

# Exact empirical AUC (Mann-Whitney), with tied scores counted as one half.
.roc_auc <- function(y_obs, y_pred) {
  n1 <- sum(y_obs == 1)
  n0 <- sum(y_obs == 0)
  r  <- rank(y_pred)
  (sum(r[y_obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Hanley-McNeil 95% confidence interval for an AUC.
.roc_hm_ci <- function(auc, n1, n0) {
  q1  <- auc / (2 - auc)
  q2  <- 2 * auc^2 / (1 + auc)
  num <- auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)
  se  <- sqrt(num / (n1 * n0))
  z   <- stats::qnorm(0.975)
  c(lo = pmax(0, auc - z * se), hi = pmin(1, auc + z * se))
}

# Everything the ROC functions report for one model: the curve to draw, the
# AUC, the Youden-optimal cut-off and the class counts.
.roc_summary <- function(y_obs, y_pred, thresholds, label) {
  if (!any(y_obs == 1) || !any(y_obs == 0))
    stop("Model '", label, "' needs both outcome classes (0 and 1) to draw ",
         "an ROC curve.", call. = FALSE)
  exact <- .roc_curve(y_obs, y_pred)
  list(
    curve = if (is.null(thresholds)) exact else
      .roc_curve(y_obs, y_pred, thresholds),
    auc   = .roc_auc(y_obs, y_pred),
    best  = exact[which.max(exact$youden), ],
    n1    = sum(y_obs == 1),
    n0    = sum(y_obs == 0)
  )
}

# Match a user-supplied outcome to a model's fitted values. Accepts an outcome
# of the same length, or one that still contains the rows the model dropped
# for missing values (removed using the model's na.action). When the model
# stores a 0/1 response, the aligned outcome must equal it.
.align_outcome <- function(outcome, model, n_pred, label) {
  na_act <- model$na.action
  if (length(outcome) != n_pred && !is.null(na_act) &&
      length(outcome) == n_pred + length(na_act))
    outcome <- outcome[-as.integer(na_act)]
  if (length(outcome) != n_pred)
    stop("Length of `outcome` (", length(outcome), ") does not match the ",
         "number of fitted values for model '", label, "' (", n_pred, ").",
         call. = FALSE)
  y_model <- model$y
  if (length(y_model) == n_pred && all(y_model %in% c(0, 1))) {
    ok <- !is.na(outcome)
    if (any(as.numeric(y_model[ok]) != outcome[ok]))
      stop("`outcome` doesn't match the response of model '", label, "'. ",
           "It must be in the same order as the data the model was ",
           "fitted to.", call. = FALSE)
  }
  outcome
}
