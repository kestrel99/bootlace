# Tests for the shared bootstrap loop (SCIENTIFIC_REVIEW.md follow-up
# "Bootstrap failures"): classify warnings, keep valid fits, record reasons.

classify_test <- function(w) {
  m <- conditionMessage(w)
  if (grepl("benign", m)) "keep" else if (grepl("separation", m)) "flag" else "drop"
}

test_that(".boot_loop() keeps, flags and drops replicates by warning class", {
  refit <- function(i) {
    switch(as.character(i %% 5),
      "0" = c(a = i, b = 1),
      "1" = { warning("benign note"); c(a = i, b = 1) },
      "2" = { warning("separation happened"); c(a = i, b = 99) },
      "3" = { warning("algorithm did not converge"); c(a = i, b = 1) },
      "4" = stop("boom")
    )
  }
  res <- .boot_loop(10L, c("a", "b"), refit, classify_test)
  expect_equal(res$n_requested, 10L)
  expect_equal(nrow(res$coefs), 6L)              # i %% 5 in 0, 1, 2
  expect_equal(res$n_flagged, 2L)
  expect_equal(sum(res$failures$n), 4L)
  expect_setequal(res$failures$reason,
                  c("warning: algorithm did not converge", "error: boom"))
})

test_that(".boot_loop() drops replicates with NA coefficients and keeps extras", {
  refit <- function(i) {
    if (i == 2L) return(list(coef = c(a = NA, b = 1), extra = "x2"))
    list(coef = c(a = i, b = 1), extra = paste0("x", i))
  }
  res <- .boot_loop(3L, c("a", "b"), refit, function(w) "drop")
  expect_equal(nrow(res$coefs), 2L)
  expect_equal(res$extra, list("x1", "x3"))
  expect_match(res$failures$reason, "NA")
})

# ---- Entry points ---------------------------------------------------------------

test_that("bootstrap_logistic() keeps refits whose only warning is non-integer weights", {
  set.seed(4)
  d   <- data.frame(age = rnorm(300, 50, 10))
  d$y <- rbinom(300, 1, plogis(-3 + 0.05 * d$age))
  d$w <- runif(300, 0.5, 2)
  fit <- suppressWarnings(glm(y ~ age, data = d, family = binomial, weights = w))
  bl  <- suppressWarnings(quiet_sample_size(
    bootstrap_logistic(fit, n_boot = 120L, seed = 1L)
  ))
  expect_equal(bl$n_boot, 120L)
  expect_equal(bl$n_requested, 120L)
})

test_that("bootstrap_logistic() keeps and counts separation replicates", {
  # Perfectly separated data: every refit warns "fitted probabilities
  # numerically 0 or 1", which used to drop every replicate.
  df_sep  <- data.frame(x = c(1:5, 20:24), y = rep(0:1, each = 5))
  fit_sep <- suppressWarnings(glm(y ~ x, data = df_sep, family = binomial))
  bl <- suppressMessages(suppressWarnings(
    bootstrap_logistic(fit_sep, n_boot = 50L, seed = 1L)
  ))
  expect_gt(bl$n_flagged, 0L)
  expect_equal(bl$n_boot + sum(bl$failures$n), bl$n_requested)
  expect_output(print(bl), "of 50 retained")
})

test_that("bootstrap_cox() keeps and counts infinite-coefficient replicates", {
  set.seed(3)
  d <- data.frame(
    time = rexp(150),
    g    = factor(c(rep("rare", 6), rep("common", 144))),
    x    = rnorm(150)
  )
  d$status <- rbinom(150, 1, 0.7)
  d$status[d$g == "rare"] <- c(1, 0, 0, 0, 0, 0)
  fit <- survival::coxph(survival::Surv(time, status) ~ x + g, data = d)
  bc  <- suppressMessages(quiet_sample_size(suppressWarnings(
    bootstrap_cox(fit, n_boot = 200L, seed = 1L)
  )))
  expect_gt(bc$n_flagged, 0L)
  expect_equal(bc$n_boot + sum(bc$failures$n), bc$n_requested)
})

test_that(".boot_loop() counts flagged replicates by labelled reason", {
  refit <- function(i) {
    if (i %% 3 == 1) warning("ridge")
    if (i %% 3 == 2) warning("hessian")
    c(a = i)
  }
  classify <- function(w) {
    m <- conditionMessage(w)
    if (m == "ridge") "flag:stopped short" else "flag:poorly identified"
  }
  res <- .boot_loop(6L, "a", refit, classify)
  expect_equal(res$n_flagged, 4L)
  expect_equal(res$flag_counts$n[match(c("stopped short", "poorly identified"),
                                       res$flag_counts$reason)], c(2L, 2L))
})
