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
# Returns c(lower, upper) with attribute "method": "BCa", or the percentile
# interval labelled with why BCa is unavailable: fewer than 3 usable
# jackknife values, or jackknife values with no spread (the acceleration is
# then undefined), or an invalid acceleration adjustment (its denominator
# 1 - a * (z0 + z) is not positive, so the adjusted probabilities are
# undefined or out of order).
.bca_ci <- function(boot, est, jack, level) {
  boot  <- boot[is.finite(boot)]
  B     <- length(boot)
  pmin_ <- 1 / (B + 1)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  # Bias correction: where the estimate sits in the bootstrap distribution.
  p0 <- (sum(boot < est) + 0.5 * sum(boot == est)) / B
  z0 <- stats::qnorm(min(max(p0, pmin_), 1 - pmin_))
  percentile <- function(why) {
    out <- stats::quantile(boot, probs, names = FALSE)
    attr(out, "method") <- paste0("percentile (BCa unavailable: ", why, ")")
    out
  }
  # Acceleration: skewness of the jackknife influence values.
  jack <- jack[is.finite(jack)]
  d    <- mean(jack) - jack
  if (length(jack) < 3L || sum(d^2) <= 0)
    return(percentile("too few usable jackknife fits"))
  a     <- sum(d^3) / (6 * sum(d^2)^1.5)
  z     <- stats::qnorm(probs)
  denom <- 1 - a * (z0 + z)
  p     <- stats::pnorm(z0 + (z0 + z) / denom)
  if (!(all(denom > 0) && all(is.finite(p)) && p[1] < p[2]))
    return(percentile("acceleration adjustment invalid"))
  out <- stats::quantile(boot, pmin(pmax(p, pmin_), 1 - pmin_), names = FALSE)
  attr(out, "method") <- "BCa"
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
#' fixed EC50 the model is a logistic GLM in the other parameters. For every
#' model, the EC50 bounds are also checked against the limiting models
#' (EC50 -> 0 and EC50 -> Inf, with any covariate effects free): if a limit is
#' within the cut-off, that bound is 0 or Inf with a warning. Other
#' profiles optimise the other parameters from two starting points and keep
#' the better fit. Before a bound is located, both ends of the bracket are
#' recomputed and checked to straddle the cut-off. Redundant (aliased)
#' covariate columns don't affect the profile. If a profile fit fails
#' numerically, that bound is `NA` with a warning, never reported as
#' infinite.
#'
#' Wald intervals are estimate +/- z x standard error from the Hessian. For
#' `ec50` they are computed on the log scale and back-transformed. They need a
#' fit with `hessian = TRUE` and a positive-definite Hessian.
#'
#' In a simulation study (100 data sets of n = 300; see
#' `tests/testthat/test-coverage.R`), nominal 95% profile intervals covered
#' E0, Emax and EC50 in 92-95% of samples in a well-identified design and
#' 89-99% in a weakly identified one, where bootstrap percentile intervals
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
        .emax_glm_profile_fn(object$x, object$y, object$Z_lin, object$Z_emax,
                             object$idx, th, nll)
      } else {
        NULL
      }
      out[k, ] <- c(
        .profile_bound(nll, th, j, -1, target, step, rescale, prof_fn),
        .profile_bound(nll, th, j, +1, target, step, rescale, prof_fn)
      )
      # The confidence set includes EC50 values arbitrarily close to 0 (or
      # Inf) whenever the limiting model is inside the cut-off, whatever the
      # walk found on the way (with or without covariates).
      if (j == object$idx$log_ec50) {
        for (side in 1:2) {
          lim <- .emax_limit_nll(object$x, object$y, object$Z_lin,
                                 object$Z_emax, object$Z_ec50, c(-1, 1)[side])
          if (!is.null(lim) && lim < target) out[k, side] <- c(-Inf, Inf)[side]
        }
      }
      for (side in which(is.na(out[k, ])))
        warning(sprintf(paste0(
          "The profile likelihood for %s could not be computed %s the ",
          "estimate (a profile fit failed numerically), so the %s bound is NA."
        ), parm[k], c("below", "above")[side], c("lower", "upper")[side]),
        call. = FALSE)
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
# it never does within `max_steps`, NA if a profile evaluation fails
# numerically (so a failure is never reported as an unbounded interval).
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
    if (is.null(p)) return(NA_real_)
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
        root <- tryCatch(
          stats::uniroot(h, ends, f.lower = fs[1], f.upper = fs[2],
                         tol = 1e-6 * max(1, abs(v)))$root,
          error = function(e) NA_real_
        )
        return(root)
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
# .profile_bound(); `starts` are unused; v may be -Inf or Inf (the limiting
# models). The value is deviance / 2, which is the negative log-likelihood
# for 0/1 outcomes. Aliased (redundant) columns get a coefficient of 0,
# which gives the same fit; NULL means the fit failed.
#
# The basis is chosen so both limits are numerically stable:
# - EC50 -> Inf: h is proportional to x, so the dose column is h scaled to a
#   maximum of 1 (a logistic fit linear in x in the limit).
# - EC50 -> 0 with no zero exposure: h = 1 - EC50 / (EC50 + x) is nearly
#   constant and numerically collinear with the intercept, which would drop
#   the dose column and give a false profile. Since h = 1 - EC50 * u with
#   u = 1 / (EC50 + x), span(1, h) = span(1, u), and u has the stable limit
#   1 / x; so below the median exposure u is used instead, and E0 and Emax
#   are recovered as e0 = a - emax, emax = -b_u / EC50.
# - EC50 -> 0 with zero exposures: h tends to the indicator of x > 0, which
#   the h column represents stably.
.emax_glm_profile_fn <- function(x, y, Z_lin, Z_emax, idx, th, nll) {
  use_u_below <- if (all(x > 0)) stats::median(x) else -Inf
  emax_pos <- seq_len(1L + length(idx$emax_cov)) + 1L
  glm_coef <- function(X) {
    f <- tryCatch(
      suppressWarnings(stats::glm.fit(X, y, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(f) || !is.finite(f$deviance)) return(NULL)
    b <- f$coefficients
    b[is.na(b)] <- 0
    list(value = f$deviance / 2, b = b)
  }
  function(v, starts) {
    ev <- exp(v)
    h  <- if (is.finite(ev)) x / (ev + x) else x / max(x)
    h[x == 0] <- 0
    if (ev < use_u_below) {
      u  <- 1 / (ev + x)
      su <- max(u)
      f  <- glm_coef(cbind(1, u / su, Z_emax * h, Z_lin))
      if (is.null(f)) return(NULL)
      b    <- f$b
      b_u  <- b[2L] / su
      emax <- if (ev > 0) -b_u / ev else -sign(b_u) * .Machine$double.xmax
      b[1L] <- b[1L] - emax
      b[2L] <- emax
    } else {
      sc <- max(h)
      if (!is.finite(sc) || sc <= 0) return(NULL)
      f <- glm_coef(cbind(1, h / sc, Z_emax * (h / sc), Z_lin))
      if (is.null(f)) return(NULL)
      b <- f$b
      # Back to the model's scale (emax and its covariate effects multiply h).
      b[emax_pos] <- b[emax_pos] / sc
      if (!is.finite(ev)) b[emax_pos] <- sign(b[emax_pos]) * .Machine$double.xmax
    }
    t <- th
    t[c(idx$e0, idx$emax0, idx$emax_cov, idx$lin)] <- b
    t[idx$log_ec50] <- v
    list(value = f$value, par = t[-idx$log_ec50])
  }
}

# Negative log-likelihood of the limiting model as the reference EC50 goes to
# 0 (`dir` = -1) or Inf (`dir` = +1), with the other parameters free:
# - EC50 -> Inf: x / (EC50_i + x) ~ x / EC50_i, so
#   eta = e0 + Z_lin b + (b0 + Z_emax g) exp(-Z_ec50 d) x;
# - EC50 -> 0, with zero exposures: h -> 1(x > 0), so
#   eta = e0 + Z_lin b + (b0 + Z_emax g) 1(x > 0);
# - EC50 -> 0, no zero exposure: h = 1 - EC50_i / (EC50_i + x), so
#   eta = a + Z_lin b + Z_emax g - (c + Z_emax[, S] e) exp(Z_ec50 d) / x,
#   where S are the Emax covariates in span(1, Z_lin) (only their effects can
#   grow without bound, offset by the linear terms).
# For fixed d each is a logistic GLM; d (the EC50 covariate effects) is
# optimised from 0. Any value found is attained in the limit, so a value
# inside a cut-off shows the confidence set reaches that boundary. NULL if
# the fit fails.
.emax_limit_nll <- function(x, y, Z_lin, Z_emax, Z_ec50, dir) {
  one <- matrix(1, length(x), 1L)
  if (dir > 0) {
    X0 <- cbind(one, Z_lin)
    W  <- cbind(one, Z_emax)
    s  <- x / max(x)
    sg <- -1
  } else if (any(x == 0)) {
    pos <- as.numeric(x > 0)
    return(.emax_glm_nll(cbind(one, Z_lin, pos, Z_emax * pos), y))
  } else {
    base   <- cbind(one, Z_lin)
    in_sp  <- vapply(seq_len(ncol(Z_emax)), function(k) {
      r <- stats::lm.fit(base, Z_emax[, k])$residuals
      max(abs(r)) <= 1e-8 * max(1, abs(Z_emax[, k]))
    }, logical(1L))
    X0 <- cbind(base, Z_emax)
    W  <- cbind(one, Z_emax[, in_sp, drop = FALSE])
    s  <- min(x) / x
    sg <- 1
  }
  dev_at <- function(d) {
    m <- if (ncol(Z_ec50) > 0L) exp(sg * drop(Z_ec50 %*% d)) else 1
    v <- .emax_glm_nll(cbind(X0, W * (s * m)), y)
    if (is.null(v)) Inf else v
  }
  if (ncol(Z_ec50) == 0L) {
    v <- dev_at(numeric(0))
    return(if (is.finite(v)) v else NULL)
  }
  o <- tryCatch(stats::optim(rep(0, ncol(Z_ec50)), dev_at, method = "BFGS"),
                error = function(e) NULL)
  v <- min(dev_at(rep(0, ncol(Z_ec50))), if (is.null(o)) Inf else o$value)
  if (is.finite(v)) v else NULL
}

# deviance / 2 of a logistic GLM (the negative log-likelihood for 0/1
# outcomes); NULL if the fit fails.
.emax_glm_nll <- function(X, y) {
  f <- tryCatch(
    suppressWarnings(stats::glm.fit(X, y, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(f) || !is.finite(f$deviance)) NULL else f$deviance / 2
}

# Scan the profile of log EC50 to classify a fit at `th` that stopped short
# of the optimum. The grid runs from 1e-4 x the lowest to 1e4 x the highest
# non-zero exposure; when its best point is at an edge, the scan continues
# outwards at the same spacing to 1e10 x, and the limiting model (EC50 -> 0
# or Inf, .emax_limit_nll(); else the far end of the scan) is computed. With `tol` = 1e-6 on the negative log-likelihood:
# - "finite": a finite EC50 (not the far end of the scan) beats `th` by more
#   than `tol`; `start` is the full parameter vector there, to restart the
#   optimiser from;
# - "diverging": the profile keeps improving to the far end of the scan and
#   the limiting model beats `th` by more than `tol`: the maximum is at the
#   EC50 boundary in direction `dir` (+1 Inf, -1 zero);
# - "flat": the best fit is towards a boundary, but nothing beats `th` by
#   more than `tol` (the fit, the best finite EC50 and the limiting model
#   are indistinguishable), or the profile turns before the far end without
#   improving on `th`: the data can't tell a very large (or small) EC50 from
#   the boundary, so EC50 is weakly identified;
# - "none": neither (an interior grid optimum no better than `th`).
.emax_boundary_check <- function(x, y, Z_lin, Z_emax, Z_ec50, idx, th, nll) {
  tol   <- 1e-6
  none  <- list(type = "none", dir = 0, start = NULL)
  x_pos <- x[x > 0]
  grid  <- seq(log(min(x_pos)) - log(1e4), log(max(x_pos)) + log(1e4),
               length.out = 41L)
  exact <- ncol(Z_ec50) == 0L
  prof <- if (exact) {
    .emax_glm_profile_fn(x, y, Z_lin, Z_emax, idx, th, nll)
  } else {
    others <- seq_along(th)[-idx$log_ec50]
    function(v, starts) {
      f <- function(r) {
        t <- th
        t[idx$log_ec50] <- v
        t[others] <- r
        nll(t)
      }
      o <- tryCatch(stats::optim(th[others], f, method = "BFGS",
                                 control = list(maxit = 1000)),
                    error = function(e) NULL)
      if (is.null(o) || !is.finite(o$value)) NULL else
        list(value = o$value, par = o$par)
    }
  }
  scan <- function(vs) {
    fits <- lapply(vs, function(v) prof(v, list()))
    list(v = vs, fits = fits,
         vals = vapply(fits, function(f) if (is.null(f)) Inf else f$value,
                       numeric(1L)))
  }
  finite_start <- function(sc, k) {
    start <- th
    start[-idx$log_ec50] <- sc$fits[[k]]$par
    start[idx$log_ec50]  <- sc$v[k]
    list(type = "finite", dir = 0, start = start)
  }
  f_th <- nll(th)
  sc   <- scan(grid)
  if (all(!is.finite(sc$vals))) return(none)
  best <- which.min(sc$vals)
  n_g  <- length(grid)
  dir  <- if (best == n_g) 1 else if (best == 1L) -1 else 0
  if (dir == 0)
    return(if (sc$vals[best] < f_th - tol) finite_start(sc, best) else none)

  # Best at an edge: continue outwards, then compare with the limit.
  spacing <- grid[2L] - grid[1L]
  edge    <- grid[if (dir > 0) n_g else 1L]
  ext     <- scan(edge + dir * spacing * seq_len(ceiling(log(1e6) / spacing)))
  sc$v    <- c(grid, ext$v)
  sc$fits <- c(sc$fits, ext$fits)
  sc$vals <- c(sc$vals, ext$vals)
  far     <- length(sc$v)
  best    <- which.min(sc$vals)
  if (best != far && sc$vals[best] < f_th - tol) return(finite_start(sc, best))
  f_lim <- .emax_limit_nll(x, y, Z_lin, Z_emax, Z_ec50, dir)
  if (is.null(f_lim)) f_lim <- sc$vals[far]
  if (best == far && f_lim < f_th - tol)
    return(list(type = "diverging", dir = dir, start = NULL))
  list(type = "flat", dir = dir, start = NULL)
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
