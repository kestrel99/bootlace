#' Bootstrap a fitted Cox proportional hazards model
#'
#' Runs case-resampling bootstrap on a fitted [survival::coxph()] object,
#' returning a matrix of bootstrapped log hazard-ratio estimates alongside the
#' original fit. The result can be inspected via [print.boot_cox()] and
#' [summary.boot_cox()].
#'
#' @details
#' Each replicate resamples the rows the original fit used and refits the
#' model from its own design: model matrix, response, case weights, offset,
#' strata and `ties` method. Rows excluded by `subset =` or by missing values
#' are never resampled.
#'
#' When the fit carries a subject or cluster identifier (`coxph(..., id =)`,
#' a `cluster =` argument, or a `cluster()` term), whole subjects are
#' resampled so that the rows of one subject stay together. Counting-process
#' (`Surv(start, stop, event)`) fits without an identifier are rejected,
#' because resampling their rows independently would break up subjects.
#'
#' The returned `data` holds the rows of the original data (from
#' `fit$call$data`) that the fit used. If the call environment cannot be
#' evaluated it falls back to `model.frame(fit)` with a warning.
#'
#' Refits that warn are not all dropped, because that would condition the
#' bootstrap on well-behaved samples. A replicate whose only warning is that a
#' coefficient may be infinite (no events in some group of the resample) is
#' kept and counted in `n_flagged`, so the percentile intervals reflect that
#' instability. Replicates that error, fail to converge, raise any other
#' warning, or have an `NA` coefficient are dropped, and the reasons are
#' counted in `failures` and reported in a message.
#'
#' @param fit A fitted [survival::coxph()] object.
#' @param n_boot Number of bootstrap replicates. Default `1000L`.
#' @param seed Integer random seed passed to [set.seed()]. Default `42L`.
#' @param conf_level Confidence level for percentile intervals, strictly
#'   between 0 and 1. Default `0.95`.
#'
#' @return An S3 object of class `"boot_cox"`, a named list with:
#'   \describe{
#'     \item{`fit`}{The original `coxph` object.}
#'     \item{`boot_coefs`}{Numeric matrix of bootstrapped log hazard-ratio
#'       estimates (`n_success` x `p`), after dropping incomplete rows.}
#'     \item{`boot_basehaz`}{List with one matrix (`time`, `hazard`) per
#'       retained replicate: that refit's uncentred baseline cumulative
#'       hazard, used by [plot_cox_survival()]. `NULL` for stratified fits.}
#'     \item{`data`}{The rows of the original data frame that the fit used.}
#'     \item{`formula`}{The model formula.}
#'     \item{`conf_level`}{The requested confidence level.}
#'     \item{`n_boot`}{Number of bootstrap replicates retained.}
#'     \item{`n_requested`}{Number of bootstrap replicates requested.}
#'     \item{`n_flagged`}{Number of retained replicates with a possibly
#'       infinite coefficient.}
#'     \item{`failures`}{Data frame (`reason`, `n`) of dropped replicates.}
#'     \item{`sample_size`}{Sample-size diagnostics; see
#'       **Sample size checks**.}
#'   }
#'
#' @seealso [summary.boot_cox()], [print.boot_cox()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L   # recode to 0/1
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex + ph.ecog,
#'     data = lung
#'   )
#'   bc <- bootstrap_cox(fit, n_boot = 200, seed = 1)
#'   print(bc)
#'   summary(bc)
#' }
#'
#' @template sample-size
#' @export
bootstrap_cox <- function(
  fit,
  n_boot     = 1000L,
  seed       = 42L,
  conf_level = 0.95
) {
  if (!inherits(fit, "coxph"))
    stop("`fit` must be a fitted coxph object.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1.")

  sample_size <- .sample_size(fit)
  .warn_sample_size(sample_size, "bootstrap_cox")

  formula <- formula(fit)

  # Refits use the fit's own design (model matrix, response, weights, offset,
  # strata, ties method), resampling whole subjects when the fit has an id.
  des <- .cox_design(fit)
  if (identical(des$type, "counting") && is.null(des$cluster))
    stop(
      "`fit` uses (start, stop] data without a subject identifier. Refit ",
      "with coxph(..., id = <subject column>) so that bootstrap_cox() can ",
      "resample whole subjects.",
      call. = FALSE
    )

  # ---- Analysis data ---------------------------------------------------------
  data <- .analysis_data(fit)
  if (is.null(data)) {
    warning(
      "bootstrap_cox: could not recover original data from fit$call$data; ",
      "falling back to model.frame(fit)."
    )
    data <- stats::model.frame(fit)
  }

  coef_names <- names(coef(fit))

  # ---- Bootstrap loop --------------------------------------------------------
  set.seed(seed)
  res <- .boot_loop(
    n_boot, coef_names,
    refit    = function(i) .refit_cox(des, .resample_rows(des), basehaz = TRUE),
    classify = .classify_cox_warning
  )
  .boot_report(res, "bootstrap_cox", "a possibly infinite coefficient")

  structure(
    list(
      fit         = fit,
      boot_coefs  = res$coefs,
      boot_basehaz = if (is.null(des$strata)) res$extra else NULL,
      data        = data,
      formula     = formula,
      conf_level  = conf_level,
      n_boot      = nrow(res$coefs),
      n_requested = res$n_requested,
      n_flagged   = res$n_flagged,
      failures    = res$failures,
      sample_size = sample_size
    ),
    class = "boot_cox"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `boot_cox` object
#'
#' Displays a compact summary of the bootstrap run and the original-data log
#' hazard-ratio estimates.
#'
#' @param x A `boot_cox` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [bootstrap_cox()], [summary.boot_cox()]
#'
#' @export
print.boot_cox <- function(x, ...) {
  cat("-- Bootstrap Cox Proportional Hazards ------------------------------\n")
  cat(sprintf("  Formula   : %s\n", deparse(x$formula)))
  cat(sprintf("  N obs     : %d\n", x$fit$n))
  if (!is.null(x$sample_size)) {
    cat(sprintf("  Events    : %s\n", .events_line(x$sample_size)))
  } else {
    cat(sprintf("  N events  : %d\n", x$fit$nevent))
  }
  cat(sprintf("  Ties      : %s\n", x$fit$method))
  cat(sprintf("  Replicates: %s\n",
              .boot_replicates_line(x, "a possibly infinite coefficient")))
  cat(sprintf("  CI level  : %.0f%%\n", x$conf_level * 100))
  cat("\nOriginal-data log hazard ratios:\n")
  print(coef(x$fit))
  invisible(x)
}


#' Summarise a `boot_cox` object
#'
#' Returns and prints a data frame of original-data estimates alongside
#' bootstrapped standard errors, bias, and percentile confidence intervals.
#' Results are presented on both the log hazard-ratio scale (for inference)
#' and the hazard-ratio scale (for interpretation).
#'
#' @param object A `boot_cox` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return A data frame (invisibly) with columns `log_HR`, `HR`, `Boot_SE`,
#'   `Boot_Bias`, `CI_<level>%_lo`, `CI_<level>%_hi`, `HR_<level>%_lo`, and
#'   `HR_<level>%_hi`, where `<level>` is `round(object$conf_level * 100)`.
#'
#' @seealso [bootstrap_cox()], [print.boot_cox()]
#'
#' @export
summary.boot_cox <- function(object, ...) {
  alpha      <- 1 - object$conf_level
  ci_pct     <- round(object$conf_level * 100)
  coef_orig  <- coef(object$fit)

  ci_lo    <- apply(object$boot_coefs, 2, stats::quantile,
                    probs = alpha / 2,       na.rm = TRUE)
  ci_hi    <- apply(object$boot_coefs, 2, stats::quantile,
                    probs = 1 - alpha / 2,   na.rm = TRUE)
  boot_se   <- apply(object$boot_coefs, 2, stats::sd,   na.rm = TRUE)
  boot_bias <- colMeans(object$boot_coefs, na.rm = TRUE) - coef_orig

  tbl <- data.frame(
    log_HR    = round(coef_orig,        4),
    HR        = round(exp(coef_orig),   4),
    Boot_SE   = round(boot_se,          4),
    Boot_Bias = round(boot_bias,        4),
    CI_lo     = round(ci_lo,            4),
    CI_hi     = round(ci_hi,            4),
    HR_lo     = round(exp(ci_lo),       4),
    HR_hi     = round(exp(ci_hi),       4),
    check.names = FALSE
  )
  names(tbl)[5:8] <- c(
    sprintf("CI_%d%%_lo",    ci_pct),
    sprintf("CI_%d%%_hi",    ci_pct),
    sprintf("HR_%d%%_lo",    ci_pct),
    sprintf("HR_%d%%_hi",    ci_pct)
  )

  cat("-- Bootstrap Cox Summary -------------------------------------------\n")
  cat(sprintf(
    "  Replicates : %d  |  CI level: %.0f%%\n\n",
    object$n_boot, object$conf_level * 100
  ))
  print(tbl)
  invisible(tbl)
}
