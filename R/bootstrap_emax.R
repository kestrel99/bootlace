#' Bootstrap a fitted logistic Emax model
#'
#' Runs case-resampling bootstrap on a fitted `logistic_emax` object, returning
#' a matrix of bootstrapped parameter vectors alongside the original fit. The
#' result can be inspected via [print.boot_emax()] and [summary.boot_emax()].
#'
#' Each replicate resamples the rows of the original data (including all
#' covariate matrices) with replacement and refits [fit_logistic_emax()] using
#' the original estimates as warm-start values. Replicates that error, whose
#' optimiser doesn't converge, or that raise any other unexpected warning are
#' dropped, and the reasons are counted in `failures` and reported in a
#' message.
#'
#' Sample-size, identifiability, diverging-ridge and EC50-range warnings from
#' the refits are not treated as failures, so replicates in which EC50 runs
#' to extreme values are kept (and counted in `flag_counts`): dropping them
#' would make the intervals look more precise than the data allow. A refit
#' whose optimiser stopped short of the optimum even after a retry, without
#' being on a diverging ridge, is an optimisation failure and is dropped;
#' see **Unrealistic EC50** in [fit_logistic_emax()]. The grouped-jackknife
#' refits used for BCa are judged the same way, and failed ones are left out
#' of the acceleration estimate. When any retained replicate has
#' EC50 outside the range the exposures can support (see **Unrealistic
#' EC50** in [fit_logistic_emax()]), one detailed warning of class
#' `bootlace_emax_ec50_warning` reports how many, the bootstrap EC50
#' distribution, and whether the percentile interval bounds are set by those
#' replicates. The sample-size check runs once, in [fit_logistic_emax()], and
#' its result is shown by [print.boot_emax()].
#'
#' @param fit A fitted [fit_logistic_emax()] object.
#' @param n_boot Number of bootstrap replicates. Default `1000L`.
#' @param seed Integer random seed passed to [set.seed()]. Default `42L`.
#' @param conf_level Confidence level for percentile intervals, strictly
#'   between 0 and 1. Default `0.95`.
#'
#' @return An S3 object of class `"boot_emax"`, a named list with:
#'   \describe{
#'     \item{`fit`}{The original `logistic_emax` object.}
#'     \item{`boot_theta`}{Numeric matrix of bootstrapped parameter vectors
#'       (`n_success` x `n_par`) in the optimisation parameterisation
#'       `(e0, emax, log_ec50, ...)`  -- the same scale as `fit$theta_hat`.}
#'     \item{`conf_level`}{The requested confidence level.}
#'     \item{`n_boot`}{Number of bootstrap replicates retained.}
#'     \item{`n_requested`}{Number of bootstrap replicates requested.}
#'     \item{`n_flagged`}{Number of retained replicates whose refit stopped
#'       short of the optimum, had a non-positive-definite Hessian, or put
#'       EC50 outside the supported range.}
#'     \item{`flag_counts`}{Data frame (`reason`, `n`): retained replicates
#'       counted by each of those reasons.}
#'     \item{`failures`}{Data frame (`reason`, `n`) of dropped replicates.}
#'   }
#'
#' @seealso [fit_logistic_emax()], [summary.boot_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#'
#' fit <- fit_logistic_emax(y, dose)
#' be  <- bootstrap_emax(fit, n_boot = 200, seed = 1)
#' print(be)
#' summary(be)
#'
#' @export
bootstrap_emax <- function(
  fit,
  n_boot     = 1000L,
  seed       = 42L,
  conf_level = 0.95
) {
  if (!inherits(fit, "logistic_emax"))
    stop("`fit` must be a fitted logistic_emax object.")
  n_boot <- .check_count(n_boot, "n_boot")
  .check_seed(seed)
  .check_prob(conf_level, "conf_level")

  n         <- length(fit$x)
  n_par     <- fit$n_par
  th_names  <- names(fit$theta_hat)

  # Warm-start vector for the jackknife refits (natural-scale ec50, as
  # fit_logistic_emax() expects). Bootstrap refits don't use it: a warm start
  # kept them near the original estimate and understated the spread of
  # EC50 by about 20% in weakly identified designs.
  warm <- as.numeric(fit$theta_hat)
  warm[fit$idx$log_ec50] <- exp(warm[fit$idx$log_ec50])

  # Pre-extract covariate matrices; NULL when 0-column (no covariates).
  Z_lin  <- if (ncol(fit$Z_lin)  > 0) fit$Z_lin  else NULL
  Z_emax <- if (ncol(fit$Z_emax) > 0) fit$Z_emax else NULL
  Z_ec50 <- if (ncol(fit$Z_ec50) > 0) fit$Z_ec50 else NULL

  refit_rows <- function(rows, start) {
    fit_logistic_emax(
      y           = fit$y[rows],
      x           = fit$x[rows],
      linear_covs = if (!is.null(Z_lin))  Z_lin[rows, , drop = FALSE]  else NULL,
      emax_covs   = if (!is.null(Z_emax)) Z_emax[rows, , drop = FALSE] else NULL,
      ec50_covs   = if (!is.null(Z_ec50)) Z_ec50[rows, , drop = FALSE] else NULL,
      start       = start,
      hessian     = FALSE,
      control     = list(maxit = 500)
    )$theta_hat
  }
  # Bootstrap refits use the automatic multi-start (start = NULL).
  refit <- function(i) refit_rows(sample.int(n, n, replace = TRUE), NULL)
  # Sample-size warnings describe the data: keep. These keep the replicate
  # but are counted under their own reason (EC50 range is also summarised
  # after the loop):
  # - EC50 diverging along the Emax/EC50 ridge, checked by the fitter: the
  #   maximum is at the EC50 -> Inf boundary. Kept on purpose: dropping them
  #   removes exactly the replicates that show EC50 is poorly determined.
  # - Hessian not positive definite (poorly identified);
  # - EC50 outside the range the exposures support.
  # Dropped, with the reason recorded: an optimiser that stopped short of the
  # optimum even after a retry and isn't on a diverging ridge (an optimisation
  # failure, so its estimates can't be trusted), optim() non-convergence, and
  # any other warning.
  classify <- function(w) {
    if (inherits(w, "bootlace_sample_size_warning")) return("keep")
    if (inherits(w, "bootlace_emax_ridge_warning"))
      return("flag:EC50 diverging along the Emax/EC50 ridge")
    if (inherits(w, "bootlace_emax_identification_warning"))
      return("flag:Hessian not positive definite (poorly identified)")
    if (inherits(w, "bootlace_emax_ec50_warning"))
      return("flag:EC50 outside the supported range")
    "drop"
  }

  set.seed(seed)
  res <- .boot_loop(n_boot, th_names, refit, classify)
  .boot_report(res, "bootstrap_emax", "weakly identified Emax/EC50")

  .warn_boot_ec50(exp(res$coefs[, fit$idx$log_ec50]), fit$x, conf_level)

  # Grouped jackknife (at most 200 groups) for the BCa acceleration used by
  # summary() and plot_emax_curve(). Leaving out a few rows moves the optimum
  # little, so these refits start from the original estimate.
  n_grp <- min(n, 200L)
  grp   <- rep_len(seq_len(n_grp), n)
  # Jackknife refits are judged like bootstrap refits: one that errors or
  # raises a "drop" warning is discarded (NA) rather than used.
  jack_theta <- t(vapply(seq_len(n_grp), function(g) {
    bad <- FALSE
    th_g <- tryCatch(
      withCallingHandlers(
        refit_rows(which(grp != g), warm),
        warning = function(w) {
          if (classify(w) == "drop") bad <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    if (is.null(th_g) || bad) rep(NA_real_, n_par) else th_g
  }, numeric(n_par)))
  colnames(jack_theta) <- th_names

  structure(
    list(
      fit         = fit,
      boot_theta  = res$coefs,
      jack_theta  = jack_theta,
      conf_level  = conf_level,
      n_boot      = nrow(res$coefs),
      n_requested = res$n_requested,
      n_flagged   = res$n_flagged,
      flag_counts = res$flag_counts,
      failures    = res$failures
    ),
    class = "boot_emax"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `boot_emax` object
#'
#' @param x A `boot_emax` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [bootstrap_emax()], [summary.boot_emax()]
#'
#' @export
print.boot_emax <- function(x, ...) {
  cat("-- Bootstrap Logistic Emax -----------------------------------------\n")
  cat(sprintf("  N obs      : %d\n",  length(x$fit$x)))
  cat(sprintf("  Parameters : %d\n",  x$fit$n_par))
  if (!is.null(x$fit$sample_size))
    cat(sprintf("  Events     : %s\n", .events_line(x$fit$sample_size)))
  cat(sprintf("  Replicates : %s\n",
              .boot_replicates_line(x, "weakly identified Emax/EC50")))
  cat(sprintf("  CI level   : %.0f%%\n", x$conf_level * 100))
  cat("\nOriginal-data estimates:\n")
  print(x$fit$coefficients)
  invisible(x)
}


#' Summarise a `boot_emax` object
#'
#' Returns and prints a data frame of original-data estimates alongside
#' bootstrapped standard errors, bias, and confidence intervals.
#'
#' By default the intervals are bias-corrected and accelerated (BCa; Efron &
#' Tibshirani, 1993). In a simulation study (100 data sets of n = 300; see
#' `tests/testthat/test-coverage.R`) the Emax estimate was biased upward and
#' skewed, and nominal 95% percentile intervals covered the true Emax in only
#' 86% of samples in a well-identified design (EC50 within the dose range)
#' and 83% in a weakly identified one (EC50 near the top dose). BCa intervals
#' covered E0, Emax and log EC50 in 90-95% and 89-93% of samples
#' respectively. Use `ci = "percentile"` for the plain percentile intervals.
#' BCa needs the jackknife fits stored by [bootstrap_emax()]; for older
#' objects without them, percentile intervals are used and the header says
#' so. For a parameter where BCa can't be computed (fewer than 3 usable
#' jackknife fits, jackknife values with no spread, or an invalid
#' acceleration adjustment) the percentile interval is shown and the output
#' says why; if that applies to every parameter, the header says percentile.
#' If no bootstrap replicates were retained, the standard errors, bias and
#' intervals are `NA` and the header says the interval is unavailable. Profile-likelihood intervals (`confint(fit)`, see
#' [confint.logistic_emax()]) did as well or better and also show when a
#' bound is not determined by the data.
#'
#' Intervals are computed on the optimisation scale (log EC50) and
#' back-transformed, so the `ec50` interval is on the natural scale.
#'
#' @param object A `boot_emax` object.
#' @param ci `"bca"` (default) or `"percentile"`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return A data frame (invisibly) with columns `Estimate`, `Boot_SE`,
#'   `Boot_Bias`, `CI_<level>%_lo`, and `CI_<level>%_hi`, and attribute
#'   `ci_method`.
#'
#' @references
#' Efron B, Tibshirani RJ (1993). *An Introduction to the Bootstrap*. Chapman
#' & Hall, New York. Chapter 14.
#'
#' @seealso [bootstrap_emax()], [print.boot_emax()], [confint.logistic_emax()]
#'
#' @export
summary.boot_emax <- function(object, ci = c("bca", "percentile"), ...) {
  ci     <- match.arg(ci)
  level  <- object$conf_level
  alpha  <- 1 - level
  ci_pct <- round(level * 100)
  fit    <- object$fit
  B      <- object$boot_theta
  th     <- fit$theta_hat
  if (ci == "bca" && is.null(object$jack_theta)) ci <- "percentile_fallback"
  # No retained replicates: no SE, bias or interval can be computed.
  if (nrow(B) == 0L) ci <- "none"

  # Intervals on the optimisation scale, then back-transform ec50.
  ints <- lapply(seq_len(ncol(B)), function(j) {
    if (ci == "none") c(NA_real_, NA_real_)
    else if (ci == "bca") .bca_ci(B[, j], th[[j]], object$jack_theta[, j], level)
    else stats::quantile(B[, j], c(alpha / 2, 1 - alpha / 2), names = FALSE,
                         na.rm = TRUE)
  })
  lims <- vapply(ints, function(r) as.numeric(r), numeric(2L))
  j_ec50 <- fit$idx$log_ec50
  lims[, j_ec50] <- exp(lims[, j_ec50])
  # Parameters where BCa was unavailable and the percentile interval is shown.
  no_bca <- if (ci == "bca") {
    names(fit$coefficients)[vapply(ints, function(r) {
      !identical(attr(r, "method"), "BCa")
    }, logical(1L))]
  } else {
    character(0)
  }

  # SE and bias on the reporting scale.
  boot_rep <- B
  boot_rep[, j_ec50] <- exp(boot_rep[, j_ec50])
  coef_orig <- fit$coefficients
  boot_se   <- apply(boot_rep, 2, stats::sd, na.rm = TRUE)
  boot_bias <- colMeans(boot_rep, na.rm = TRUE) - coef_orig
  boot_bias[!is.finite(boot_bias)] <- NA_real_

  tbl <- data.frame(
    Estimate  = round(coef_orig, 4),
    Boot_SE   = round(boot_se,   4),
    Boot_Bias = round(boot_bias, 4),
    CI_lo     = round(lims[1, ], 4),
    CI_hi     = round(lims[2, ], 4),
    row.names = names(coef_orig),
    check.names = FALSE
  )
  names(tbl)[4:5] <- c(
    sprintf("CI_%d%%_lo", ci_pct),
    sprintf("CI_%d%%_hi", ci_pct)
  )
  method_lbl <- switch(ci,
    none = sprintf("unavailable (no bootstrap replicates retained, 0 of %d requested)",
                   max(0L, object$n_requested)),
    bca = "BCa",
    percentile = "percentile",
    percentile_fallback = "percentile (no jackknife fits stored for BCa)"
  )
  # BCa requested but unavailable for every parameter: say so in the header.
  if (ci == "bca" && length(no_bca) == length(coef_orig)) {
    reasons <- unique(vapply(ints, function(r) {
      sub("^percentile \\(BCa unavailable: (.*)\\)$", "\\1", attr(r, "method"))
    }, character(1L)))
    method_lbl <- sprintf("percentile (BCa unavailable: %s)",
                          paste(reasons, collapse = "; "))
  }
  attr(tbl, "ci_method")      <- method_lbl
  attr(tbl, "bca_unavailable") <- no_bca

  cat("-- Bootstrap Emax Summary ------------------------------------------\n")
  cat(sprintf(
    "  Replicates : %d  |  CI level: %.0f%%  |  CI method: %s\n",
    object$n_boot, level * 100, method_lbl
  ))
  if (length(no_bca) > 0L) {
    why <- vapply(ints, function(r) {
      sub("^percentile \\(BCa unavailable: (.*)\\)$", "\\1", attr(r, "method"))
    }, character(1L))[match(no_bca, names(fit$coefficients))]
    for (w in unique(why))
      cat(sprintf(
        "  BCa unavailable (%s) for: %s; percentile intervals shown.\n",
        w, paste(no_bca[why == w], collapse = ", ")
      ))
  }
  cat("\n")
  print(tbl)
  invisible(tbl)
}


# Warn, in detail, when bootstrap EC50 estimates fall outside the observed
# non-zero exposure range, and say whether the percentile interval is
# affected.
.warn_boot_ec50 <- function(ec50, x, conf_level) {
  lim   <- .ec50_support(x)
  hi    <- lim[["hi"]]
  lo    <- lim[["lo"]]
  above <- sum(ec50 > hi)
  below <- sum(ec50 < lo)
  n     <- length(ec50)
  if (n == 0L || above + below == 0L) return(invisible(NULL))

  tail_share <- (1 - conf_level) / 2
  ci_hit <- c(
    if (above / n > tail_share) "upper",
    if (below / n > tail_share) "lower"
  )
  q <- stats::quantile(ec50, c(0.5, tail_share, 1 - tail_share), names = FALSE)
  .emax_warning(paste0(
    sprintf(
      "bootstrap_emax(): EC50 lies outside the range the exposures can support (%s, a tenth of the lowest non-zero exposure, to %s, the highest) in %d of %d bootstrap replicates (%.1f%%): %d above and %d below. ",
      format(lo), format(hi), above + below, n, 100 * (above + below) / n,
      above, below
    ),
    sprintf(
      "Bootstrap EC50: median %s, %.0f%% percentile interval %s to %s, largest %s. ",
      format(signif(q[1], 3)), 100 * conf_level, format(signif(q[2], 3)),
      format(signif(q[3], 3)), format(signif(max(ec50), 3))
    ),
    "In those replicates the data don't pin down where the response levels ",
    "off, so Emax and EC50 trade off against each other and take extreme ",
    "values together. ",
    if (length(ci_hit) > 0L) {
      sprintf(paste0(
        "This is more than %.1f%% of replicates, so the %s bound(s) of the ",
        "percentile intervals for EC50 and Emax are set by these unsupported ",
        "values and shouldn't be reported as estimates. "
      ), 100 * tail_share, paste(ci_hit, collapse = " and "))
    } else {
      sprintf(paste0(
        "This is no more than %.1f%% of replicates in either direction, so ",
        "the %.0f%% percentile intervals are not set by these values, but ",
        "EC50 is only weakly determined. "
      ), 100 * tail_share, 100 * conf_level)
    },
    "Options: a profile-likelihood interval for EC50 (confint(fit, method = ",
    "\"profile\")), fixing Emax or EC50 ",
    "from prior information, a simpler exposure-response model, or data at ",
    "more informative exposures."
  ), "bootlace_emax_ec50_warning")
}
