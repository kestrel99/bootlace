#' Bootstrap a fitted binomial GLM
#'
#' Runs case-resampling bootstrap on a fitted `glm` object with
#' `family = binomial`, returning a matrix of bootstrapped coefficients
#' alongside the original fit. The result is consumed by
#' [plot_logistic_curves()] to construct pointwise confidence bands, and
#' can also be inspected directly via [print.boot_logistic()] and
#' [summary.boot_logistic()].
#'
#' @details
#' Each replicate resamples the rows the original fit used and refits the
#' model from its own design: model matrix, response, prior weights and
#' offset. Rows excluded by `subset =` or by missing values are never
#' resampled, and case weights and offsets are kept in every refit.
#'
#' For plotting, the function also recovers the untransformed rows of the
#' original data (from `fit$call$data`, restricted to the rows the fit used),
#' so that inline formula transformations (e.g. `log(age)`) can be applied to
#' prediction grids. If the call environment cannot be evaluated --for example
#' when the fit was constructed inside another function --it falls back to
#' `model.frame(fit)` with a warning; in that case formulae containing inline
#' transformations may not plot correctly.
#'
#' Bootstrap replicates that error or generate warnings during refitting are
#' skipped and counted. If any replicates are skipped, a message is emitted.
#' Replicates with any `NA` coefficient (e.g. arising from complete separation
#' in a resample) are also dropped before the result is returned.
#'
#' @param fit A fitted [glm()] object with `family = binomial` (any link).
#' @param n_boot Number of bootstrap replicates. Default `1000L`.
#' @param seed Integer random seed passed to [set.seed()] for reproducibility.
#'   Default `42L`.
#' @param conf_level Confidence level for percentile intervals, strictly
#'   between 0 and 1. Default `0.95`.
#'
#' @return An S3 object of class `"boot_logistic"`, a named list with elements:
#'   \describe{
#'     \item{`fit`}{The original `glm` object.}
#'     \item{`boot_coefs`}{Numeric matrix of bootstrapped coefficients from
#'       successful refits (`n_success` x `p`), after removing skipped
#'       replicates and rows with any `NA` coefficient.}
#'     \item{`data`}{The rows of the original (untransformed) data frame that
#'       the fit used.}
#'     \item{`formula`}{The model formula.}
#'     \item{`conf_level`}{The requested confidence level.}
#'     \item{`n_boot`}{Number of successful bootstrap replicates retained.}
#'     \item{`sample_size`}{Sample-size diagnostics; see
#'       **Sample size checks**.}
#'   }
#'
#' @seealso [plot_logistic_curves()], [summary.boot_logistic()]
#'
#' @examples
#' df <- data.frame(
#'   y   = rbinom(200, 1, 0.4),
#'   age = rnorm(200, 50, 10),
#'   sex = factor(sample(c("M", "F"), 200, replace = TRUE))
#' )
#' fit <- glm(y ~ age + sex, data = df, family = binomial)
#' bl  <- bootstrap_logistic(fit, n_boot = 200, seed = 1)
#' print(bl)
#' summary(bl)
#'
#' @template sample-size
#' @export
bootstrap_logistic <- function(
  fit,
  n_boot = 1000L,
  seed = 42L,
  conf_level = 0.95
) {
  # ---- Input validation ------------------------------------------------------
  if (!inherits(fit, "glm")) {
    stop("`fit` must be a fitted glm object.")
  }
  if (fit$family$family != "binomial") {
    stop("`fit` must use family = binomial (any link).")
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be strictly between 0 and 1.")
  }

  sample_size <- .sample_size(fit)
  .warn_sample_size(sample_size, "bootstrap_logistic")

  # ---- Analysis data ---------------------------------------------------------
  # The untransformed rows the fit used, so the plotting functions can apply
  # inline formula transformations (e.g. log(age)) to their own grids.
  formula <- formula(fit)
  data <- .analysis_data(fit)
  if (is.null(data)) {
    warning(
      "bootstrap_logistic: could not recover original data from ",
      "fit$call$data; falling back to model.frame(fit). ",
      "Inline formula transformations (e.g. log(x)) may not plot correctly."
    )
    data <- model.frame(fit)
  }

  # Refits use the fit's own design (model matrix, response, prior weights,
  # offset), so the analysis rows, weights and offsets can't change.
  des <- .glm_design(fit)
  coef_names <- names(coef(fit))

  # ---- Bootstrap loop --------------------------------------------------------
  set.seed(seed)
  boot_coefs <- matrix(
    NA_real_,
    nrow = n_boot,
    ncol = length(coef_names),
    dimnames = list(NULL, coef_names)
  )
  n_failed <- 0L

  for (i in seq_len(n_boot)) {
    idx <- .resample_rows(des)

    warned <- FALSE
    b <- tryCatch(
      withCallingHandlers(
        .refit_glm(des, idx),
        warning = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )

    if (is.null(b) || warned) {
      n_failed <- n_failed + 1L
      next
    }

    # Aliased columns (e.g. a factor level absent from the resample) are NA
    boot_coefs[i, ] <- b
  }

  if (n_failed > 0L) {
    message(sprintf(
      "bootstrap_logistic: %d replicate(s) were skipped due to errors or warnings.",
      n_failed
    ))
  }

  # Drop replicates with any NA coefficient (e.g. from complete separation)
  boot_coefs <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]

  n_success <- nrow(boot_coefs)
  if (n_success < 100L) {
    warning(sprintf(
      "Only %d bootstrap replicates succeeded; CIs may be unreliable.",
      n_success
    ))
  }

  structure(
    list(
      fit = fit,
      boot_coefs = boot_coefs,
      data = data,
      formula = formula,
      conf_level = conf_level,
      n_boot = n_success,
      sample_size = sample_size
    ),
    class = "boot_logistic"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `boot_logistic` object
#'
#' Displays a compact summary of the bootstrap run and the original-data
#' coefficients.
#'
#' @param x A `boot_logistic` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [bootstrap_logistic()], [summary.boot_logistic()]
#'
#' @export
print.boot_logistic <- function(x, ...) {
  cat("-- Bootstrap Logistic Regression ----------------------------------\n")
  cat(sprintf("  Formula   : %s\n", deparse(x$formula)))
  cat(sprintf("  Link      : %s\n", x$fit$family$link))
  cat(sprintf("  N obs     : %d\n", nrow(x$data)))
  if (!is.null(x$sample_size))
    cat(sprintf("  Events    : %s\n", .events_line(x$sample_size)))
  cat(sprintf("  Replicates: %d successful\n", x$n_boot))
  cat(sprintf("  CI level  : %.0f%%\n", x$conf_level * 100))
  cat("\nOriginal-data coefficients:\n")
  print(coef(x$fit))
  invisible(x)
}


#' Summarise a `boot_logistic` object
#'
#' Returns and prints a data frame of original-data estimates alongside
#' bootstrapped standard errors, bias, and percentile confidence intervals
#' for each model coefficient.
#'
#' @param object A `boot_logistic` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return A data frame (invisibly) with columns `Estimate`, `Boot_SE`,
#'   `Boot_Bias`, `CI_<level>%_lo`, and `CI_<level>%_hi`, where `<level>` is
#'   `round(object$conf_level * 100)`.
#'
#' @seealso [bootstrap_logistic()], [print.boot_logistic()]
#'
#' @export
summary.boot_logistic <- function(object, ...) {
  alpha <- 1 - object$conf_level
  coef_orig <- coef(object$fit)
  ci_lo <- apply(
    object$boot_coefs,
    2,
    stats::quantile,
    probs = alpha / 2,
    na.rm = TRUE
  )
  ci_hi <- apply(
    object$boot_coefs,
    2,
    stats::quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  boot_se <- apply(object$boot_coefs, 2, stats::sd, na.rm = TRUE)
  boot_bias <- colMeans(object$boot_coefs, na.rm = TRUE) - coef_orig

  ci_pct <- round(object$conf_level * 100)
  tbl <- data.frame(
    Estimate = round(coef_orig, 4),
    Boot_SE = round(boot_se, 4),
    Boot_Bias = round(boot_bias, 4),
    CI_lo = round(ci_lo, 4),
    CI_hi = round(ci_hi, 4),
    check.names = FALSE
  )
  names(tbl)[4:5] <- c(
    sprintf("CI_%d%%_lo", ci_pct),
    sprintf("CI_%d%%_hi", ci_pct)
  )

  cat("-- Bootstrap Summary -----------------------------------------------\n")
  cat(sprintf(
    "  Replicates : %d  |  CI level: %.0f%%\n\n",
    object$n_boot,
    object$conf_level * 100
  ))
  print(tbl)
  invisible(tbl)
}
