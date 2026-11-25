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
