# Emax interval methods: BCa bootstrap intervals, multi-start bootstrap
# refits, and profile-likelihood / Wald confint(). Motivated by the coverage
# study (test-coverage.R), where percentile intervals under-covered.

emax_set <- function(ec50 = 10, emax = 3, n = 300L, seed = 1L) {
  set.seed(seed)
  x <- rep(c(0, 5, 10, 25, 50, 100), length.out = n)
  list(x = x, y = rbinom(n, 1, plogis(qlogis(0.1) + emax * x / (ec50 + x))))
}
quiet_emax <- function(expr) suppressMessages(muffle_class(
  expr, c("bootlace_emax_ec50_warning", "bootlace_emax_identification_warning")
))

# ---- BCa ------------------------------------------------------------------------

test_that(".bca_ci() is the percentile interval when bias and acceleration are 0", {
  boot <- seq(-1, 1, length.out = 1000)        # symmetric, none equal to 0
  expect_equal(.bca_ci(boot, 0, c(-1, 1), 0.9),
               unname(quantile(boot, c(0.05, 0.95))))
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
