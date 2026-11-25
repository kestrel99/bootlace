# Emax interval methods: BCa bootstrap intervals, multi-start bootstrap
# refits, and profile-likelihood / Wald confint(). Motivated by the coverage
# study (test-coverage.R), where percentile intervals under-covered.

emax_set <- function(ec50 = 10, emax = 3, n = 300L, seed = 1L) {
  set.seed(seed)
  x <- rep(c(0, 5, 10, 25, 50, 100), length.out = n)
  list(x = x, y = rbinom(n, 1, plogis(qlogis(0.1) + emax * x / (ec50 + x))))
}
quiet_emax <- function(expr) suppressMessages(muffle_class(
  expr, c("bootlace_emax_ec50_warning", "bootlace_emax_identification_warning",
          "bootlace_emax_stopped_warning")
))

# ---- BCa ------------------------------------------------------------------------

test_that(".bca_ci() is the percentile interval when bias and acceleration are 0", {
  boot <- seq(-1, 1, length.out = 1000)        # symmetric, none equal to 0
  expect_equal(.bca_ci(boot, 0, c(-1, 1), 0.9),
               unname(quantile(boot, c(0.05, 0.95))), ignore_attr = TRUE)
})

test_that(".bca_ci() moves the interval towards a biased estimate", {
  set.seed(1)
  boot <- rnorm(4000, mean = 1)                # bootstrap centred above est
  ci   <- .bca_ci(boot, 0, c(-1, 1), 0.95)
  expect_true(all(ci < quantile(boot, c(0.025, 0.975))))
})

test_that("bootstrap_emax() stores jackknife fits; summary() defaults to BCa", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  be  <- quiet_emax(bootstrap_emax(fit, n_boot = 120L, seed = 1L))
  expect_equal(dim(be$jack_theta), c(200L, fit$n_par))   # min(n, 200) groups

  out <- capture.output(tb <- summary(be))
  expect_match(out, "BCa", all = FALSE)

  capture.output(tp <- summary(be, ci = "percentile"))
  expect_equal(tp[["CI_95%_lo"]][2],
               round(unname(quantile(be$boot_theta[, "emax"], 0.025)), 4))
  expect_false(isTRUE(all.equal(tb, tp)))
})

test_that("bootstrap_emax() refits use automatic multi-start", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  starts   <- logical(0)
  real_fit <- fit_logistic_emax
  local_mocked_bindings(fit_logistic_emax = function(..., start = NULL) {
    starts <<- c(starts, is.null(start))
    real_fit(..., start = start)
  })
  quiet_emax(bootstrap_emax(fit, n_boot = 110L, seed = 1L))
  # bootstrap refits: automatic starts; jackknife refits: warm start
  expect_equal(sum(starts), 110L)
})

# ---- confint() --------------------------------------------------------------------

test_that("profile intervals contain the estimate and agree with Wald when well identified", {
  d   <- emax_set(n = 600L)
  fit <- fit_logistic_emax(d$y, d$x)
  pr  <- confint(fit, method = "profile")
  wd  <- confint(fit, method = "wald")
  expect_equal(rownames(pr), names(fit$coefficients))
  expect_true(all(pr[, 1] < fit$coefficients & fit$coefficients < pr[, 2]))
  expect_equal(pr["e0", ], wd["e0", ], tolerance = 0.1)
})

test_that("Wald intervals come from the Hessian, with EC50 on the log scale", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  wd  <- confint(fit, method = "wald", level = 0.9)
  se  <- sqrt(diag(fit$vcov_theta))
  z   <- qnorm(0.95)
  expect_equal(unname(wd["e0", ]),
               unname(fit$theta_hat[["e0"]] + c(-1, 1) * z * se[["e0"]]))
  expect_equal(unname(wd["ec50", ]),
               unname(exp(fit$theta_hat[["log_ec50"]] + c(-1, 1) * z * se[["log_ec50"]])))
})

test_that("an unidentified EC50 has an infinite profile upper bound and a warning", {
  d   <- emax_set(ec50 = 1000, emax = 4, n = 400L, seed = 2L)
  fit <- quiet_emax(fit_logistic_emax(d$y, d$x))
  expect_warning(ci <- confint(fit, parm = "ec50", method = "profile"),
                 "upper")
  expect_equal(unname(ci["ec50", 2]), Inf)
})

test_that("EC50 warnings point to the profile interval", {
  d   <- emax_set(ec50 = 1000, emax = 4, n = 400L, seed = 2L)
  msg <- NULL
  withCallingHandlers(
    fit_logistic_emax(d$y, d$x),
    warning = function(w) {
      if (inherits(w, "bootlace_emax_ec50_warning")) msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_match(msg, "confint(fit, method = \"profile\")", fixed = TRUE)
})

# ---- Curve ribbon -----------------------------------------------------------------

test_that("plot_emax_curve() offers BCa (default) and percentile ribbons", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  be  <- quiet_emax(bootstrap_emax(fit, n_boot = 120L, seed = 1L))
  pb  <- plot_emax_curve(be, n_grid = 20L)
  pp  <- plot_emax_curve(be, n_grid = 20L, ci = "percentile")
  b   <- ggplot2::ggplot_build(pb)
  rib_b <- Filter(function(z) "ymin" %in% names(z), b$data)[[1]]
  rib_p <- Filter(function(z) "ymin" %in% names(z), ggplot2::ggplot_build(pp)$data)[[1]]
  expect_false(isTRUE(all.equal(rib_b$ymin, rib_p$ymin)))
  expect_match(pb$labels$subtitle, "BCa")
})

# ---- Follow-up review: optimiser failures vs poor identification -----------------

test_that("stopping short and a non-PD Hessian are different warning classes", {
  d <- emax_set()
  first_class <- function() {
    cls <- character(0)
    withCallingHandlers(
      fit_logistic_emax(d$y, d$x),
      warning = function(w) {
        cls <<- c(cls, class(w)[1])
        invokeRestart("muffleWarning")
      }
    )
    cls
  }
  # Optimiser stopped short: PD Hessian but a large Newton decrement.
  local_mocked_bindings(.polish_optimum = function(f, p, ...) {
    list(par = p, H_inv = diag(length(p)), decrement = 1)
  })
  expect_true("bootlace_emax_stopped_warning" %in% first_class())
  # Not identified: no positive-definite Hessian.
  local_mocked_bindings(.polish_optimum = function(f, p, ...) {
    list(par = p, H_inv = NULL, decrement = NA_real_)
  })
  expect_true("bootlace_emax_identification_warning" %in% first_class())
})

test_that("bootstrap_emax() keeps ridge replicates but counts each reason", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  be  <- quiet_emax(bootstrap_emax(fit, n_boot = 120L, seed = 1L))
  expect_named(be$flag_counts, c("reason", "n"))
  expect_equal(sum(grepl("stopped short|identified|EC50", be$flag_counts$reason)),
               nrow(be$flag_counts))
  expect_equal(be$n_boot + sum(be$failures$n), be$n_requested)
})

# ---- Follow-up review round 2 -------------------------------------------------------

# Independent profile for EC50 (no covariates): at fixed EC50 the model is a
# logistic GLM in x / (EC50 + x).
glm_profile_ci <- function(fit, level = 0.95) {
  x <- fit$x
  y <- fit$y
  prof <- function(le) {
    h <- x / (exp(le) + x)
    deviance(suppressWarnings(glm(y ~ h, family = binomial))) / 2
  }
  target <- -fit$loglik + qchisq(level, 1) / 2
  le_hat <- log(fit$coefficients[["ec50"]])
  grid   <- le_hat + seq(-15, 15, by = 0.05)
  above  <- vapply(grid, prof, numeric(1)) > target
  lo_i   <- max(c(0, which(above & grid < le_hat)))
  hi_i   <- min(c(Inf, which(above & grid > le_hat)))
  c(if (lo_i == 0) 0 else exp(grid[lo_i]), if (is.infinite(hi_i)) Inf else exp(grid[hi_i]))
}

test_that("confint() doesn't fail on converged fits with awkward profiles", {
  for (s in c(83L, 113L)) {
    set.seed(s)
    x <- rep(c(0, 5, 10, 25, 50, 100), length.out = 200)
    y <- rbinom(200, 1, plogis(qlogis(0.1) + 3 * x / (100 + x)))
    fit <- quiet_emax(fit_logistic_emax(y, x))
    ci  <- suppressWarnings(confint(fit, parm = "ec50"))
    ref <- glm_profile_ci(fit)
    # same bound type (finite / infinite), and finite bounds agree to the
    # 0.05 log-scale grid of the reference
    expect_equal(unname(is.finite(ci[1, ]) & ci[1, ] > 0),
                 is.finite(ref) & ref > 0)
    ok <- is.finite(ref) & ref > 0
    expect_equal(log(unname(ci[1, ok])), log(ref[ok]), tolerance = 0.06)
  }
})

test_that("BCa falls back to a labelled percentile interval past the pole", {
  set.seed(1)
  boot <- rexp(1000)
  # Acceleration about +0.16 (one low jackknife outlier) and the estimate far
  # in the bootstrap tail (z0 about 2): 1 - a * (z0 + z) < 0 at level 0.99999.
  jack <- c(rep(0, 99), -60)
  ci   <- .bca_ci(boot, unname(quantile(boot, 0.98)), jack, 0.99999)
  expect_true(ci[1] <= ci[2])
  expect_match(attr(ci, "method"), "percentile")
  ok <- .bca_ci(boot, 1, c(-1, 1), 0.95)
  expect_identical(attr(ok, "method"), "BCa")
})

test_that("a refit that still stops short is told apart from a diverging ridge", {
  d <- emax_set()
  classes <- function() {
    cls <- character(0)
    withCallingHandlers(
      fit_logistic_emax(d$y, d$x),
      warning = function(w) {
        cls <<- c(cls, class(w)[1])
        invokeRestart("muffleWarning")
      }
    )
    cls
  }
  local_mocked_bindings(.polish_optimum = function(f, p, ...) {
    list(par = p, H_inv = diag(length(p)), decrement = 1)
  })
  local_mocked_bindings(.emax_ridge_diverges = function(...) TRUE)
  expect_true("bootlace_emax_ridge_warning" %in% classes())
  local_mocked_bindings(.emax_ridge_diverges = function(...) FALSE)
  expect_true("bootlace_emax_stopped_warning" %in% classes())
})

test_that("bootstrap_emax() drops stopped-short refits and keeps ridge ones", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  real_fit <- fit_logistic_emax
  calls <- 0L
  local_mocked_bindings(fit_logistic_emax = function(..., start = NULL) {
    calls <<- calls + 1L
    if (is.null(start) && calls %% 4L == 0L)
      .emax_warning("stopped", "bootlace_emax_stopped_warning")
    if (is.null(start) && calls %% 4L == 1L)
      .emax_warning("ridge", "bootlace_emax_ridge_warning")
    real_fit(..., start = start)
  })
  be <- suppressMessages(muffle_class(
    bootstrap_emax(fit, n_boot = 140L, seed = 1L), "bootlace_emax_ec50_warning"
  ))
  expect_equal(sum(be$failures$n[grepl("stopped", be$failures$reason)]), 35L)
  expect_true(any(grepl("ridge", be$flag_counts$reason)))
})

test_that("jackknife refits that stop short are discarded, not used", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  real_fit <- fit_logistic_emax
  local_mocked_bindings(fit_logistic_emax = function(..., start = NULL) {
    if (!is.null(start)) .emax_warning("stopped", "bootlace_emax_stopped_warning")
    real_fit(..., start = start)
  })
  be <- quiet_emax(bootstrap_emax(fit, n_boot = 110L, seed = 1L))
  expect_true(all(is.na(be$jack_theta)))
})

test_that("summary() says when BCa was unavailable and percentile is shown", {
  d   <- emax_set()
  fit <- fit_logistic_emax(d$y, d$x)
  be  <- quiet_emax(bootstrap_emax(fit, n_boot = 110L, seed = 1L))
  local_mocked_bindings(.bca_ci = function(boot, est, jack, level) {
    structure(quantile(boot, c(0.025, 0.975), names = FALSE),
              method = "percentile (BCa unavailable)")
  })
  out <- capture.output(tb <- summary(be))
  expect_match(out, "BCa unavailable", all = FALSE)
  expect_equal(attr(tb, "bca_unavailable"), names(fit$coefficients))
})
