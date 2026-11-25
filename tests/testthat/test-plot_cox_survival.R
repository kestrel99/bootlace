# Tests for plot_cox_survival(), including SCIENTIFIC_REVIEW.md follow-up
# "Cox survival bands": the baseline hazard is re-estimated per replicate.

lung_boot <- function(n_boot = 120L) {
  lung <- survival::lung
  lung$status <- lung$status - 1L
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + sex, data = lung
  )
  suppressMessages(bootstrap_cox(fit, n_boot = n_boot, seed = 1L))
}

test_that("bootstrap_cox() stores one baseline hazard per retained replicate", {
  bc <- lung_boot()
  expect_length(bc$boot_basehaz, bc$n_boot)
  expect_equal(colnames(bc$boot_basehaz[[1]]), c("time", "hazard"))
})

test_that("survival bands use each replicate's own baseline hazard", {
  bc  <- lung_boot()
  nd  <- data.frame(age = 60, sex = 1)
  t0  <- 365
  p   <- plot_cox_survival(bc, nd, times = t0)
  row <- p$data[p$data$time == t0, ]

  h0  <- vapply(bc$boot_basehaz, function(h) {
    stats::approx(h[, "time"], h[, "hazard"], xout = t0, method = "constant",
                  f = 0, yleft = 0, rule = 2)$y
  }, numeric(1))
  lp  <- drop(as.matrix(nd[, c("age", "sex")]) %*% t(bc$boot_coefs))
  s_b <- exp(-h0 * exp(lp))
  alpha <- 1 - bc$conf_level
  expect_equal(row$ci_lo, unname(quantile(s_b, alpha / 2)))
  expect_equal(row$ci_hi, unname(quantile(s_b, 1 - alpha / 2)))
})

test_that("survival is 1 before the first event time", {
  bc <- lung_boot()
  p  <- plot_cox_survival(bc, data.frame(age = 60, sex = 1), times = c(0, 1))
  expect_true(all(p$data$surv == 1))
  expect_true(all(p$data$ci_lo == 1 & p$data$ci_hi == 1))
})

test_that("objects without per-replicate baselines fall back with a note", {
  bc <- lung_boot()
  bc$boot_basehaz <- NULL
  p <- plot_cox_survival(bc, data.frame(age = 60, sex = 1))
  expect_match(p$labels$subtitle, "fixed baseline")
})

test_that("interval bounds are ordered on every scale (follow-up review blocker)", {
  bc <- lung_boot()
  nd <- data.frame(age = 60, sex = 1)
  for (f in c("survival", "event", "cumhaz")) {
    d <- plot_cox_survival(bc, nd, times = c(180, 365, 730), fun = f)$data
    expect_true(all(d$ci_lo <= d$surv & d$surv <= d$ci_hi), label = f)
  }
  # cumhaz bounds are -log of the survival bounds, in swapped order (up to
  # quantile() interpolation, which doesn't commute exactly with -log)
  s <- plot_cox_survival(bc, nd, times = 365)$data
  h <- plot_cox_survival(bc, nd, times = 365, fun = "cumhaz")$data
  expect_equal(h$ci_lo, -log(s$ci_hi), tolerance = 1e-4)
  expect_equal(h$ci_hi, -log(s$ci_lo), tolerance = 1e-4)
})
