# Shared bootstrap loop for bootstrap_logistic(), bootstrap_cox() and
# bootstrap_emax().
#
# Dropping every replicate that warns would condition the bootstrap on
# well-behaved samples. Instead each warning is classified:
#   "keep" - benign; the replicate is kept
#   "flag" - informative (e.g. separation); the replicate is kept and counted
#   "drop" - the fit can't be trusted (e.g. non-convergence); the replicate is
#            dropped and the reason recorded
# Errors and NA coefficients also drop the replicate, with a reason.

# `refit(i)` returns a named coefficient vector, or list(coef =, extra =)
# when the caller needs more per replicate (kept in `extra` for retained
# replicates). `classify(w)` returns "keep", "flag" or "drop".
.boot_loop <- function(n_boot, coef_names, refit, classify) {
  coefs   <- matrix(NA_real_, n_boot, length(coef_names),
                    dimnames = list(NULL, coef_names))
  extra   <- vector("list", n_boot)
  reason  <- rep(NA_character_, n_boot)
  flagged <- logical(n_boot)

  for (i in seq_len(n_boot)) {
    drop_msg <- NULL
    flag     <- FALSE
    res <- tryCatch(
      withCallingHandlers(
        refit(i),
        warning = function(w) {
          cls <- classify(w)
          if (cls == "drop" && is.null(drop_msg))
            drop_msg <<- paste("warning:", conditionMessage(w))
          if (cls == "flag") flag <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        drop_msg <<- paste("error:", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(drop_msg)) {
      reason[i] <- drop_msg
      next
    }
    b <- if (is.list(res)) res$coef else res
    if (anyNA(b)) {
      reason[i] <- "NA coefficient (e.g. a factor level absent from the resample)"
      next
    }
    coefs[i, ] <- b
    if (is.list(res)) extra[i] <- list(res$extra)
    flagged[i] <- flag
  }

  keep <- is.na(reason)
  failures <- if (any(!keep)) {
    tab <- sort(table(reason[!keep]), decreasing = TRUE)
    data.frame(reason = names(tab), n = as.integer(tab),
               stringsAsFactors = FALSE)
  } else {
    data.frame(reason = character(0), n = integer(0), stringsAsFactors = FALSE)
  }
  list(
    coefs       = coefs[keep, , drop = FALSE],
    extra       = extra[keep],
    n_requested = as.integer(n_boot),
    n_flagged   = sum(flagged[keep]),
    failures    = failures
  )
}

# Message listing dropped replicates by reason, and the number kept despite
# an informative warning (`flag_what` describes it).
.boot_report <- function(res, caller, flag_what) {
  n_drop <- sum(res$failures$n)
  if (n_drop > 0L)
    message(sprintf(
      "%s: %d of %d replicate(s) dropped: %s.",
      caller, n_drop, res$n_requested,
      paste0(res$failures$reason, " (", res$failures$n, ")", collapse = "; ")
    ))
  if (res$n_flagged > 0L)
    message(sprintf(
      "%s: %d retained replicate(s) showed %s; they are kept so the percentile intervals reflect that instability.",
      caller, res$n_flagged, flag_what
    ))
  if (nrow(res$coefs) < 100L)
    warning(sprintf(
      "Only %d bootstrap replicates succeeded; CIs may be unreliable.",
      nrow(res$coefs)
    ), call. = FALSE)
  invisible(NULL)
}

# Print line, e.g. "980 of 1000 retained (12 flagged: separation)".
.boot_replicates_line <- function(x, flag_what) {
  if (is.null(x$n_requested)) return(sprintf("%d successful", x$n_boot))
  out <- sprintf("%d of %d retained", x$n_boot, x$n_requested)
  if (isTRUE(x$n_flagged > 0L))
    out <- sprintf("%s (%d with %s)", out, x$n_flagged, flag_what)
  out
}

# Warning classifiers ------------------------------------------------------------

.classify_glm_warning <- function(w) {
  m <- conditionMessage(w)
  if (grepl("non-integer #successes", m, fixed = TRUE)) return("keep")
  if (grepl("fitted probabilities numerically 0 or 1", m, fixed = TRUE))
    return("flag")
  "drop"
}

.classify_cox_warning <- function(w) {
  if (grepl("coefficient may be infinite", conditionMessage(w), fixed = TRUE))
    return("flag")
  "drop"
}
