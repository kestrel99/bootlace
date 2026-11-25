# Argument validation (SCIENTIFIC_REVIEW.md follow-up "Argument validation").

lung01 <- function() {
  d <- survival::lung
  d$status <- d$status - 1L
  d
}

test_that("bootstrap functions reject bad counts, seeds and confidence levels", {
  for (bad in list(0, -5, 2.5, NA, "10", c(10, 20))) {
    expect_error(bootstrap_logistic(FIT, n_boot = bad), "`n_boot`")
  }
  expect_error(bootstrap_logistic(FIT, seed = NA), "`seed`")
  expect_error(bootstrap_logistic(FIT, conf_level = NA), "`conf_level`")
  expect_error(bootstrap_logistic(FIT, conf_level = "0.9"), "`conf_level`")
  fit_c <- survival::coxph(survival::Surv(time, status) ~ age, data = lung01())
  expect_error(bootstrap_cox(fit_c, n_boot = 0), "`n_boot`")
  set.seed(1)
  x <- rep(c(0, 5, 10, 25, 50, 100), length.out = 200)
  y <- rbinom(200, 1, plogis(-2 + 3 * x / (10 + x)))
  fit_e <- fit_logistic_emax(y, x)
  expect_error(bootstrap_emax(fit_e, n_boot = 1.5), "`n_boot`")
  expect_error(bootstrap_emax(fit_e, conf_level = 1), "`conf_level`")
})

test_that("vpc_cox() validates its simulation arguments", {
  d   <- lung01()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  expect_error(vpc_cox(fit, d, n_sim = 0), "`n_sim`")
  expect_error(vpc_cox(fit, d, n_grid = 1), "`n_grid`")
  expect_error(vpc_cox(fit, d, pi_probs = c(0.95, 0.05)), "`pi_probs`")
  expect_error(vpc_cox(fit, d, pi_probs = c(0, 0.95)), "`pi_probs`")
  expect_error(vpc_cox(fit, d, pi_probs = 0.5), "`pi_probs`")
  expect_error(vpc_cox(fit, d, t_max = -1), "`t_max`")
  expect_error(vpc_cox(fit, d, seed = "a"), "`seed`")
})

test_that("plotting and ROC functions validate times, thresholds and grids", {
  d   <- lung01()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  bc  <- suppressMessages(bootstrap_cox(fit, n_boot = 100L, seed = 1L))
  nd  <- data.frame(age = 60, sex = 1)
  expect_error(plot_cox_survival(bc, nd, times = c(10, Inf)), "`times`")
  expect_error(plot_cox_survival(bc, nd, times = -1), "`times`")
  expect_error(roc_logistic(list(m = FIT), DF$y, thresholds = 1.5), "`thresholds`")
  expect_error(plot_logistic_curves(BL, n_grid = 1), "`n_grid`")
})

test_that("predict.logistic_emax() validates exposure and level", {
  set.seed(1)
  x <- rep(c(0, 5, 10, 25, 50, 100), length.out = 200)
  y <- rbinom(200, 1, plogis(-2 + 3 * x / (10 + x)))
  fit <- fit_logistic_emax(y, x)
  expect_error(predict(fit, newdata = c(1, -2)), "non-negative")
  expect_error(predict(fit, newdata = c(1, Inf)), "finite")
  expect_error(predict(fit, newdata = 10, type = "ci", level = 2), "`level`")
})

test_that("stepwise_reduce() and fp_data() validate their arguments", {
  expect_error(stepwise_reduce(FIT, alpha = NA), "`alpha`")
  expect_error(stepwise_reduce(FIT, max_iterations = 0), "`max_iterations`")
  expect_error(fp_data(FIT, DF, ci_level = 0), "`ci_level`")
  expect_error(fp_data(FIT, DF, percs = c(0.5, 1.2)), "`percs`")
})

test_that("compare_logistic_roc() rejects NA and infinite n_cuts / auc_digits", {
  m2 <- glm(y ~ age, data = DF, family = binomial)
  for (bad in list(NA, Inf, -1, 1.5, "2")) {
    expect_error(compare_logistic_roc(list(FIT, m2), n_cuts = bad), "`n_cuts`")
    expect_error(compare_logistic_roc(list(FIT, m2), auc_digits = bad),
                 "`auc_digits`")
  }
})
