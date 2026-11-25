# Tests for fit_logistic_emax(), including SCIENTIFIC_REVIEW.md finding 3:
# no convergence claims or precise standard errors for unidentified models.

emax_data <- function(n = 300L, seed = 1L, doses = c(0, 5, 10, 25, 50, 100),
                      e0 = qlogis(0.1), emax = 3, ec50 = 10) {
  set.seed(seed)
  x <- rep(doses, length.out = n)
  list(x = x, y = rbinom(n, 1, plogis(e0 + emax * x / (ec50 + x))))
}

# ---- Input validation -----------------------------------------------------------

test_that("all-zero exposure is rejected instead of reported as converged", {
  d <- emax_data()
  expect_error(fit_logistic_emax(d$y, rep(0, length(d$y))), "distinct")
})

test_that("fewer than 3 distinct exposure values are rejected", {
  d <- emax_data(doses = c(0, 10))
  expect_error(fit_logistic_emax(d$y, d$x), "distinct")
})

test_that("negative or infinite exposure is rejected", {
  d <- emax_data()
  x <- d$x
  x[1] <- -1
  expect_error(fit_logistic_emax(d$y, x), "non-negative")
  x[1] <- Inf
  expect_error(fit_logistic_emax(d$y, x), "finite")
})

test_that("a single outcome class is rejected", {
  d <- emax_data()
  expect_error(fit_logistic_emax(rep(1, length(d$x)), d$x), "both outcome")
})

test_that("invalid starting values are rejected", {
  d <- emax_data()
  expect_error(fit_logistic_emax(d$y, d$x, start = c(0, 1, -5)), "start")
  expect_error(fit_logistic_emax(d$y, d$x, start = c(NA, 1, 5)), "start")
})

# ---- Covariance -------------------------------------------------------------------

test_that(".pd_inverse() refuses singular, indefinite or non-finite matrices", {
  expect_null(.pd_inverse(matrix(1, 2, 2)))
  expect_null(.pd_inverse(diag(c(1, -1))))
  expect_null(.pd_inverse(matrix(c(1, NA, NA, 1), 2)))
  expect_null(.pd_inverse(diag(c(1, 1e-14))))
  expect_equal(.pd_inverse(diag(c(2, 4))), diag(c(0.5, 0.25)))
})

# ---- Estimation ---------------------------------------------------------------------

test_that("a well-identified design converges with finite standard errors", {
  d   <- emax_data()
  fit <- fit_logistic_emax(d$y, d$x)
  expect_true(fit$converged)
  expect_false(is.null(fit$vcov_theta))
  expect_true(all(is.finite(fit$std_err) & fit$std_err > 0))
})

test_that("parameters are recovered from a large simulated trial", {
  d   <- emax_data(n = 6000L, seed = 7L)
  fit <- fit_logistic_emax(d$y, d$x)
  est <- fit$theta_hat[c("e0", "emax", "log_ec50")]
  se  <- sqrt(diag(fit$vcov_theta))[c("e0", "emax", "log_ec50")]
  truth <- c(qlogis(0.1), 3, log(10))
  expect_true(all(abs(est - truth) < 4 * se))
})

test_that("multiple starts don't do worse than the default start", {
  d     <- emax_data(seed = 3L)
  auto  <- fit_logistic_emax(d$y, d$x)
  fixed <- fit_logistic_emax(d$y, d$x, start = c(qlogis(0.2), 1, 25))
  expect_gte(auto$loglik, fixed$loglik - 1e-6)
})
