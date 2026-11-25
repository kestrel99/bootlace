# Regression tests for SCIENTIFIC_REVIEW.md finding 2: ROC outcome alignment
# and AUC calculation.

make_roc_df <- function(n = 200L, seed = 3L) {
  set.seed(seed)
  d <- data.frame(age = rnorm(n, 55, 10), grp = factor(sample(1:3, n, TRUE)))
  d$y <- rbinom(n, 1, plogis(-4 + 0.07 * d$age))
  d
}

# Exact empirical AUC: P(score of a case > score of a control) + 0.5 * P(tie).
exact_auc <- function(y, p) {
  pos <- p[y == 1]
  neg <- p[y == 0]
  mean(outer(pos, neg, ">") + 0.5 * outer(pos, neg, "=="))
}

test_that("roc_logistic() AUC equals the exact empirical AUC", {
  d   <- make_roc_df()
  fit <- glm(y ~ age, data = d, family = binomial)
  roc <- roc_logistic(list(m = fit), outcome = d$y)
  expect_equal(unname(roc$auc), exact_auc(d$y, fitted(fit)), tolerance = 1e-12)
})

test_that("roc_logistic() AUC accounts for tied scores", {
  d   <- make_roc_df()
  fit <- glm(y ~ grp, data = d, family = binomial)   # 3 distinct scores
  roc <- roc_logistic(list(m = fit), outcome = d$y)
  expect_equal(unname(roc$auc), exact_auc(d$y, fitted(fit)), tolerance = 1e-12)
})

test_that("plotting thresholds don't change the AUC or the optimal cut-off", {
  d    <- make_roc_df()
  fit  <- glm(y ~ age, data = d, family = binomial)
  full <- roc_logistic(list(m = fit), outcome = d$y)
  one  <- roc_logistic(list(m = fit), outcome = d$y, thresholds = 0.5)
  expect_equal(one$auc, full$auc)
  expect_equal(one$best_threshold, full$best_threshold)
  expect_true(full$best_threshold$m$threshold %in% fitted(fit))
})

test_that("roc_logistic() rejects an outcome of the wrong length", {
  d   <- make_roc_df(n = 7L)
  fit <- suppressWarnings(glm(y ~ age, data = d, family = binomial))
  expect_error(roc_logistic(list(m = fit), outcome = c(d$y, 0, 1)), "Length")
})

test_that("roc_logistic() aligns an outcome that still has the dropped rows", {
  d <- make_roc_df()
  d$age[c(20, 90, 150)] <- NA
  fit <- glm(y ~ age, data = d, family = binomial)
  roc <- roc_logistic(list(m = fit), outcome = d$y)
  expect_equal(unname(roc$auc), exact_auc(fit$y, fitted(fit)), tolerance = 1e-12)
})

test_that("roc_logistic() rejects an outcome in a different order", {
  d   <- make_roc_df()
  fit <- glm(y ~ age, data = d, family = binomial)
  expect_error(roc_logistic(list(m = fit), outcome = rev(d$y)), "order")
})

test_that("roc_logistic() needs both outcome classes", {
  d <- make_roc_df()
  d$y[] <- 1
  fit <- suppressWarnings(glm(y ~ age, data = d, family = binomial))
  expect_error(roc_logistic(list(m = fit), outcome = d$y), "both")
})

test_that("roc_logistic_emax() AUC is exact and independent of thresholds", {
  set.seed(42)
  dose <- rep(c(0, 5, 10, 25, 50, 100), length.out = 200)
  y    <- rbinom(200, 1, plogis(qlogis(0.1) + 3 * dose / (10 + dose)))
  fit  <- fit_logistic_emax(y, dose)
  full <- roc_logistic_emax(fit)
  one  <- roc_logistic_emax(fit, thresholds = 0.5)
  expect_equal(unname(full$auc), exact_auc(fit$y, fit$fitted.values),
               tolerance = 1e-12)
  expect_equal(one$auc, full$auc)
})
