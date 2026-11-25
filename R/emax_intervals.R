# Interval methods for logistic Emax models: the shared negative
# log-likelihood, BCa bootstrap intervals, and confint() with profile-
# likelihood or Wald intervals.

# Negative log-likelihood of the logistic Emax model as a function of theta
# (e0, emax, log_ec50, linear..., emax covariates..., ec50 covariates...),
# with `idx` as stored in a logistic_emax object.
.emax_nll_fn <- function(x, y, Z_lin, Z_emax, Z_ec50, idx) {
  p_lin  <- ncol(Z_lin)
  p_emax <- ncol(Z_emax)
  p_ec50 <- ncol(Z_ec50)
  function(theta) {
    emax_i <- theta[idx$emax0] +
      if (p_emax > 0) drop(Z_emax %*% theta[idx$emax_cov]) else 0
    ec50_i <- exp(theta[idx$log_ec50] +
      if (p_ec50 > 0) drop(Z_ec50 %*% theta[idx$ec50_cov]) else 0)
    eta <- theta[idx$e0] + emax_i * x / (ec50_i + x) +
      if (p_lin > 0) drop(Z_lin %*% theta[idx$lin]) else 0
    -sum(
      y       * stats::plogis(eta, log.p = TRUE) +
      (1 - y) * stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
    )
  }
}

# Bias-corrected and accelerated (BCa) bootstrap interval for one statistic:
# `boot` are its bootstrap replicates, `est` its estimate on the original
# data, `jack` its (grouped) jackknife estimates. Efron & Tibshirani (1993),
# An Introduction to the Bootstrap, ch. 14.
#
# Returns c(lower, upper) with attribute "method": "BCa", or, when the
# acceleration adjustment is invalid (its denominator 1 - a * (z0 + z) is not
# positive, so the adjusted probabilities are undefined or out of order),
# the percentile interval labelled "percentile (BCa unavailable)".
.bca_ci <- function(boot, est, jack, level) {
  boot  <- boot[is.finite(boot)]
  B     <- length(boot)
  pmin_ <- 1 / (B + 1)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  # Bias correction: where the estimate sits in the bootstrap distribution.
  p0 <- (sum(boot < est) + 0.5 * sum(boot == est)) / B
  z0 <- stats::qnorm(min(max(p0, pmin_), 1 - pmin_))
  # Acceleration: skewness of the jackknife influence values.
  jack <- jack[is.finite(jack)]
  d    <- mean(jack) - jack
  a    <- if (length(jack) > 2L && sum(d^2) > 0) sum(d^3) / (6 * sum(d^2)^1.5) else 0
  z     <- stats::qnorm(probs)
  denom <- 1 - a * (z0 + z)
  p     <- stats::pnorm(z0 + (z0 + z) / denom)
  ok    <- all(denom > 0) && all(is.finite(p)) && p[1] < p[2]
  out <- if (ok) {
    stats::quantile(boot, pmin(pmax(p, pmin_), 1 - pmin_), names = FALSE)
  } else {
    stats::quantile(boot, probs, names = FALSE)
  }
  attr(out, "method") <- if (ok) "BCa" else "percentile (BCa unavailable)"
  out
}


#' Confidence intervals for a logistic Emax model
#'
#' Profile-likelihood (default) or Wald confidence intervals for the
#' parameters of a [fit_logistic_emax()] fit.
#'
#' Profile-likelihood intervals collect the parameter values whose profile
#' log-likelihood (maximised over the other parameters) is within
#' `qchisq(level, 1) / 2` of the maximum. Unlike Wald intervals they don't
#' assume a symmetric, quadratic likelihood, which matters for EC50 and Emax
#' when the exposures don't pin down the plateau. When the profile never
#' falls far enough in one direction, that bound is infinite (0 for the lower
#' bound of `ec50`), and a warning says so: the data don't determine the
#' parameter in that direction.
#'
#' For EC50 in a model without EC50 covariates the profile is exact: at a
#' fixed EC50 the model is a logistic GLM in the other parameters. Other
#' profiles optimise the other parameters from two starting points and keep
#' the better fit. Before a bound is located, both ends of the bracket are
#' recomputed and checked to straddle the cut-off.
#'
#' Wald intervals are estimate +/- z x standard error from the Hessian. For
#' `ec50` they are computed on the log scale and back-transformed. They need a
#' fit with `hessian = TRUE` and a positive-definite Hessian.
#'
#' In a simulation study (100 data sets of n = 300; see
#' `tests/testthat/test-coverage.R`), nominal 95% profile intervals covered
#' E0, Emax and EC50 in 92-95% of samples in a well-identified design and
#' 91-99% in a weakly identified one, where bootstrap percentile intervals
#' covered Emax in only 86% and 83%. The BCa intervals of
#' [summary.boot_emax()] covered 89-95%.
#'
#' @param object A `logistic_emax` object.
#' @param parm Parameters to report, by name (as in `coef(object)`) or
#'   position. Default: all.
#' @param level Confidence level. Default `0.95`.
#' @param method `"profile"` (default) or `"wald"`.
#' @param ... Currently unused.
#'
#' @return A matrix with one row per parameter and columns for the lower and
#'   upper bounds, on the reporting scale of `coef(object)` (`ec50` on the
#'   natural scale, EC50 covariate effects on the log scale).
#'
#' @seealso [fit_logistic_emax()], [summary.boot_emax()]
#'
#' @examples
#' set.seed(42)
#' dose <- rep(c(0, 5, 10, 25, 50, 100), length.out = 300)
#' y    <- rbinom(300, 1, plogis(qlogis(0.1) + 3 * dose / (10 + dose)))
#' fit  <- fit_logistic_emax(y, dose)
#' confint(fit)
#' confint(fit, method = "wald")
#'
#' @export
confint.logistic_emax <- function(object, parm, level = 0.95,
                                  method = c("profile", "wald"), ...) {
  method <- match.arg(method)
  .check_prob(level, "level")
  coef_nm <- names(object$coefficients)
  if (missing(parm)) parm <- coef_nm
  if (is.numeric(parm)) parm <- coef_nm[parm]
  bad <- setdiff(parm, coef_nm)
  if (length(bad) > 0L)
    stop("Unknown parameter(s) in `parm`: ", paste(bad, collapse = ", "),
         call. = FALSE)

  # Coefficients and theta share positions; only ec50 is on another scale.
  j_all <- match(parm, coef_nm)
  th    <- object$theta_hat
  a     <- (1 - level) / 2
  out   <- matrix(NA_real_, length(parm), 2L, dimnames = list(
    parm, paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  ))

  se <- if (is.null(object$vcov_theta)) NULL else sqrt(diag(object$vcov_theta))
  if (method == "wald") {
    if (is.null(se))
      stop("Wald intervals need a positive-definite Hessian; this fit has ",
           "none (refit with hessian = TRUE, or use method = \"profile\").",
           call. = FALSE)
    z <- stats::qnorm(1 - a)
    for (k in seq_along(j_all)) {
      j <- j_all[k]
      out[k, ] <- th[[j]] + c(-1, 1) * z * se[[j]]
    }
  } else {
    nll    <- .emax_nll_fn(object$x, object$y, object$Z_lin, object$Z_emax,
                           object$Z_ec50, object$idx)
    target <- nll(th) + stats::qchisq(level, 1) / 2
    for (k in seq_along(j_all)) {
      j    <- j_all[k]
      step <- if (!is.null(se) && is.finite(se[[j]]) && se[[j]] > 0) se[[j]] else 0.25
      # Along the Emax/EC50 ridge the two grow together; keep the other one's
      # starting value on the ridge as this one moves.
      rescale <- if (j %in% c(object$idx$log_ec50, object$idx$emax0)) {
        .emax_ridge_rescale(object, j)
      } else {
        NULL
      }
      # At a fixed EC50 the model is a logistic GLM in the other parameters
      # (when EC50 has no covariates), so profile log EC50 exactly.
      prof_fn <- if (j == object$idx$log_ec50 && ncol(object$Z_ec50) == 0L) {
        .emax_glm_profile(object, nll)
      } else {
        NULL
      }
      out[k, ] <- c(
        .profile_bound(nll, th, j, -1, target, step, rescale, prof_fn),
        .profile_bound(nll, th, j, +1, target, step, rescale, prof_fn)
      )
      for (side in which(is.infinite(out[k, ])))
        warning(sprintf(paste0(
          "The profile likelihood for %s doesn't fall to the %g%% cut-off %s ",
          "the estimate, so the %s bound is %s: the data don't determine %s ",
          "in that direction."
        ), parm[k], 100 * level, c("below", "above")[side],
        c("lower", "upper")[side],
        if (j == object$idx$log_ec50 && side == 1L) "0" else
          c("-Inf", "Inf")[side], parm[k]), call. = FALSE)
    }
  }

  is_ec50 <- j_all == object$idx$log_ec50
  out[is_ec50, ] <- exp(out[is_ec50, ])
  out
}

# Where the profile negative log-likelihood of parameter `j` reaches `target`,
# walking from the estimate in direction `dir` with growing steps; +/-Inf if
# it never does.
# `rescale(v, v_prev, start)`, if given, adjusts the other parameters'
# starting values when parameter `j` moves from `v_prev` to `v`.
# `profile_fn(v, starts)`, if given, computes the profile exactly; otherwise
# the other parameters are optimised from each of `starts` and the best kept.
# Before root-finding, both ends of the bracket are recomputed with the same
# function and must straddle `target`; if the far end turns out to be below
# it after all (a better nuisance fit), the walk continues from there.
.profile_bound <- function(nll, th, j, dir, target, step, rescale = NULL,
                           profile_fn = NULL, max_steps = 40L) {
  others <- seq_along(th)[-j]
  if (is.null(profile_fn)) {
    profile_fn <- function(v, starts) {
      f <- function(r) {
        t <- th
        t[j] <- v
        t[others] <- r
        nll(t)
      }
      best <- NULL
      for (st in starts) {
        o <- tryCatch(
          stats::optim(st, f, method = "BFGS",
                       control = list(maxit = 1000, reltol = 1e-12)),
          error = function(e) NULL
        )
        if (!is.null(o) && is.finite(o$value) &&
            (is.null(best) || o$value < best$value))
          best <- list(value = o$value, par = o$par)
      }
      best
    }
  }
  v_prev   <- th[[j]]
  par_prev <- th[others]
  for (k in seq_len(max_steps)) {
    v     <- v_prev + dir * step
    start <- if (is.null(rescale)) par_prev else rescale(v, v_prev, par_prev)
    p <- profile_fn(v, list(start, par_prev))
    if (is.null(p)) break
    if (p$value >= target) {
      starts <- list(par_prev, p$par)
      h <- function(u) {
        q <- profile_fn(u, starts)
        if (is.null(q)) NA_real_ else q$value - target
      }
      h_near <- h(v_prev)
      h_far  <- h(v)
      if (is.finite(h_near) && is.finite(h_far) && h_near < 0 && h_far >= 0) {
        ends <- if (dir > 0) c(v_prev, v) else c(v, v_prev)
        fs   <- if (dir > 0) c(h_near, h_far) else c(h_far, h_near)
        return(stats::uniroot(h, ends, f.lower = fs[1], f.upper = fs[2],
                              tol = 1e-6 * max(1, abs(v)))$root)
      }
      # Not a verified crossing: carry on from the better fit at v.
      if (is.finite(h_far) && h_far < 0) {
        q <- profile_fn(v, starts)
        if (!is.null(q)) p <- q
      }
    }
    v_prev   <- v
    par_prev <- p$par
    step     <- step * 1.5
  }
  dir * Inf
}

# Exact profile of log EC50 for a model without EC50 covariates: at fixed
# EC50 the linear predictor is e0 + (emax + Z_emax g) h + Z_lin b with
# h = x / (EC50 + x), a logistic GLM. Returns profile_fn(v, starts) for
# .profile_bound(); `starts` are unused.
.emax_glm_profile <- function(object, nll) {
  idx <- object$idx
  x   <- object$x
  y   <- object$y
  th  <- object$theta_hat
  function(v, starts) {
    h <- x / (exp(v) + x)
    X <- cbind(1, h, object$Z_emax * h, object$Z_lin)
    f <- tryCatch(
      suppressWarnings(stats::glm.fit(X, y, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(f) || anyNA(f$coefficients)) return(NULL)
    t <- th
    t[c(idx$e0, idx$emax0, idx$emax_cov, idx$lin)] <- f$coefficients
    t[idx$log_ec50] <- v
    list(value = nll(t), par = t[-idx$log_ec50])
  }
}

# Start-value adjuster for profiling log EC50 or Emax (parameter `j`) that
# keeps Emax / (EC50 + x_m) about constant, x_m the median non-zero exposure,
# so the curve at typical exposures barely changes between profile steps:
# - profiling log EC50: scale the Emax parameters by
#   (EC50 + x_m) / (EC50_prev + x_m);
# - profiling Emax: move EC50 so that EC50 + x_m scales with Emax.
.emax_ridge_rescale <- function(object, j) {
  x_m    <- stats::median(object$x[object$x > 0])
  idx    <- object$idx
  others <- seq_along(object$theta_hat)[-j]
  if (j == idx$log_ec50) {
    pos <- match(c(idx$emax0, idx$emax_cov), others)
    return(function(v, v_prev, start) {
      start[pos] <- start[pos] * (exp(v) + x_m) / (exp(v_prev) + x_m)
      start
    })
  }
  pos <- match(idx$log_ec50, others)
  function(v, v_prev, start) {
    if (v_prev == 0 || v / v_prev <= 0) return(start)
    ec50 <- (exp(start[pos]) + x_m) * v / v_prev - x_m
    if (ec50 > 0) start[pos] <- log(ec50)
    start
  }
}
