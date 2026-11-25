# Sample-size diagnostics: events per variable (EPV) and related checks.
#
# Every model-fitting or bootstrap entry point calls .sample_size() on the fit
# it works with, stores the result as `sample_size`, and passes it to
# .warn_sample_size(). Thresholds live here so they are changed in one place.

.EPV_MIN       <- 10  # Peduzzi et al. (1996): problems below 10
.EPV_SELECTION <- 25  # Heinze, Wallisch & Dunkler (2018), Table 3
.EVENTS_MIN    <- 20  # pragmatic floor, not a literature standard
.X_LEVELS_MIN  <- 4   # distinct exposures (incl. zero) for more than an exactly determined curve
# Significance level equivalent to AIC for a 1-df test (about 0.157). Heinze
# et al. (2018, Section 2.4.1) advise no smaller alpha when EPV < 25.
.ALPHA_AIC     <- stats::pchisq(2, df = 1, lower.tail = FALSE)

# Compute sample-size diagnostics for a binomial glm, coxph, or logistic_emax
# fit. `selection = TRUE` applies the variable-selection EPV tiers, with
# `alpha` the selection significance level.
.sample_size <- function(fit, selection = FALSE, alpha = NULL) {
  ss <- if (inherits(fit, "coxph")) {
    .counts_cox(fit)
  } else if (inherits(fit, "logistic_emax")) {
    .counts_emax(fit)
  } else {
    .counts_glm(fit)
  }
  ss$epv <- if (ss$n_params > 0) ss$limiting_events / ss$n_params else Inf
  ss$epv_threshold <- if (selection) .EPV_SELECTION else .EPV_MIN
  ss$problems <- .sample_size_problems(ss, selection, alpha)
  ss
}

.counts_glm <- function(fit) {
  w <- fit$prior.weights
  b <- stats::coef(fit)
  has_int <- "(Intercept)" %in% names(b) && !is.na(b[["(Intercept)"]])
  mf <- stats::model.frame(fit)
  list(
    type            = "binary",
    n               = length(fit$y),
    events          = sum(w * fit$y),
    non_events      = sum(w * (1 - fit$y)),
    limiting_events = min(sum(w * fit$y), sum(w * (1 - fit$y))),
    n_params        = sum(!is.na(b)) - has_int,
    sparse_levels   = .sparse_levels(.predictor_cols(mf), w * fit$y, w * (1 - fit$y)),
    n_x_levels      = NA_integer_
  )
}

.counts_cox <- function(fit) {
  y <- fit$y
  if (is.null(y)) y <- stats::model.response(stats::model.frame(fit))
  status <- y[, ncol(y)]
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  list(
    type            = "survival",
    n               = fit$n,
    events          = fit$nevent,
    non_events      = NA_real_,
    limiting_events = fit$nevent,
    n_params        = sum(!is.na(stats::coef(fit))),
    sparse_levels   = if (is.null(mf)) character(0)
                      else .sparse_levels(.predictor_cols(mf), status),
    n_x_levels      = NA_integer_
  )
}

.counts_emax <- function(fit) {
  y <- fit$y
  x_lv <- unique(fit$x[is.finite(fit$x)])
  groups <- c(
    .emax_cov_groups(fit$Z_lin,  fit$mm_info_lin,  "linear_covs"),
    .emax_cov_groups(fit$Z_emax, fit$mm_info_emax, "emax_covs"),
    .emax_cov_groups(fit$Z_ec50, fit$mm_info_ec50, "ec50_covs")
  )
  # The same covariate can enter more than one parameter; check it once.
  groups <- groups[!duplicated(names(groups))]
  list(
    type            = "binary",
    n               = length(y),
    events          = sum(y),
    non_events      = sum(1 - y),
    limiting_events = min(sum(y), sum(1 - y)),
    n_params        = fit$n_par - 1L,
    sparse_levels   = .sparse_levels(groups, y, 1 - y),
    n_x_levels      = length(x_lv)
  )
}

# Predictor columns of a model frame: drops the response, "(weights)"-style
# extras, and strata()/cluster() terms, which are not estimated parameters.
.predictor_cols <- function(mf) {
  keep <- seq_along(mf) != 1L &
    !startsWith(names(mf), "(") &
    !startsWith(names(mf), "strata(") &
    !startsWith(names(mf), "cluster(")
  mf[keep]
}

# Rebuild grouping variables from an Emax covariate matrix: one factor per
# original categorical covariate (via mm_info), plus every 0/1 numeric column
# that isn't one of those dummies.
.emax_cov_groups <- function(Z, mm_info, label) {
  if (ncol(Z) == 0L) return(list())
  out <- list()
  used <- character(0)
  for (v in names(mm_info$xlevels)) {
    lv <- mm_info$xlevels[[v]]
    cols <- paste0(v, lv[-1])
    if (!all(cols %in% colnames(Z))) next
    D <- Z[, cols, drop = FALSE]
    lvl_idx <- ifelse(rowSums(D) == 0, 1L, max.col(D, ties.method = "first") + 1L)
    name <- if (v == ".var") label else v
    out[[name]] <- factor(lv[lvl_idx], levels = lv)
    used <- c(used, cols)
  }
  for (j in setdiff(colnames(Z), used)) {
    z <- Z[, j]
    if (all(z %in% c(0, 1))) out[[j]] <- factor(z, levels = c(0, 1))
  }
  out
}

# Levels of categorical predictors with no events (or, when `non_events` is
# given, no non-events). Returns labels like "sex = F (0 events)".
.sparse_levels <- function(preds, events, non_events = NULL) {
  out <- character(0)
  for (v in names(preds)) {
    g <- preds[[v]]
    if (!(is.factor(g) || is.character(g) || is.logical(g))) next
    g <- droplevels(as.factor(g))
    ev <- tapply(events, g, sum)
    out <- c(out, sprintf("%s = %s (0 events)", v, names(ev)[ev == 0]))
    if (is.null(non_events)) next
    ne <- tapply(non_events, g, sum)
    out <- c(out, sprintf("%s = %s (0 non-events)", v, names(ne)[ne == 0]))
  }
  out
}

.sample_size_problems <- function(ss, selection = FALSE, alpha = NULL) {
  epv_txt <- sprintf(
    "EPV = %.1f (%s events, %d parameters)",
    ss$epv, format(ss$limiting_events), ss$n_params
  )
  probs <- character(0)
  # Selection tiers follow Heinze et al. (2018), Table 3: EPV <= 10 and
  # 10 < EPV <= 25.
  if (selection && ss$epv <= .EPV_MIN) {
    probs <- c(probs, paste0(
      epv_txt, " is at most ", .EPV_MIN,
      "; variable selection is not recommended. Fit the full model instead,",
      " ideally with shrinkage or penalised likelihood."
    ))
  } else if (selection && ss$epv <= .EPV_SELECTION) {
    msg <- paste0(
      epv_txt, " is at most ", .EPV_SELECTION,
      "; selected coefficients are likely biased away from zero. Follow",
      " selection with shrinkage or use penalised estimation, and check",
      " selection stability with bootstrap."
    )
    if (!is.null(alpha) && alpha < .ALPHA_AIC)
      msg <- paste0(msg, sprintf(
        " At this EPV, alpha = %g is too strict; use at least %.3f (AIC).",
        alpha, .ALPHA_AIC
      ))
    probs <- c(probs, msg)
  } else if (!selection && ss$epv < .EPV_MIN) {
    probs <- c(probs, paste0(
      epv_txt, " is below ", .EPV_MIN,
      "; estimates are likely to be biased and overfitted."
    ))
  }
  if (ss$limiting_events < .EVENTS_MIN) {
    what <- if (ss$type == "survival") "events" else
      "events in the less frequent outcome"
    probs <- c(probs, sprintf(
      "Only %s %s (fewer than %d); intervals may be unreliable even for a simple model.",
      format(ss$limiting_events), what, .EVENTS_MIN
    ))
  }
  if (length(ss$sparse_levels) > 0L) {
    probs <- c(probs, paste0(
      "Sparse factor levels: ", paste(ss$sparse_levels, collapse = ", "),
      "; separation or infinite estimates are likely."
    ))
  }
  if (!is.na(ss$n_x_levels) && ss$n_x_levels < .X_LEVELS_MIN) {
    probs <- c(probs, sprintf(
      "Only %d distinct exposure levels (counting zero); with fewer than 4, E0, Emax and EC50 are at best exactly determined and EC50 is usually poorly estimated.",
      ss$n_x_levels
    ))
  }
  probs
}

# Raise one bootlace_sample_size_warning listing every problem in `ss`.
.warn_sample_size <- function(ss, caller) {
  if (length(ss$problems) == 0L) return(invisible(ss))
  msg <- paste0(
    caller, "(): the data may not support this model.\n",
    paste0("* ", ss$problems, collapse = "\n")
  )
  cond <- structure(
    class = c("bootlace_sample_size_warning", "warning", "condition"),
    list(message = msg, call = NULL)
  )
  warning(cond)
  invisible(ss)
}

# One-line summary for print methods, e.g. "54 events, 66 non-events (EPV 10.8, 5 params)".
.events_line <- function(ss) {
  epv <- if (is.finite(ss$epv)) sprintf("%.1f", ss$epv) else "n/a"
  counts <- if (ss$type == "survival") {
    sprintf("%s of %d", format(ss$events), ss$n)
  } else {
    sprintf("%s events, %s non-events", format(ss$events), format(ss$non_events))
  }
  sprintf("%s (EPV %s, %d params)", counts, epv, ss$n_params)
}
