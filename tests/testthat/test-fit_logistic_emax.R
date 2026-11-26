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

# ---- Unrealistic EC50 and bootstrap replicates ---------------------------------

# Run `expr`, muffling every warning; return the value and the warnings.
collect_warnings <- function(expr) {
  ws  <- list()
  val <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}
has_class <- function(ws, cls) any(vapply(ws, inherits, logical(1L), cls))

test_that("an EC50 above the highest exposure gives a detailed warning", {
  d   <- emax_data(n = 400L, seed = 2L, emax = 4, ec50 = 1000)
  res <- collect_warnings(fit_logistic_emax(d$y, d$x))
  expect_true(has_class(res$warnings, "bootlace_emax_ec50_warning"))
  msg <- conditionMessage(Filter(
    function(w) inherits(w, "bootlace_emax_ec50_warning"), res$warnings
  )[[1]])
  expect_match(msg, "above the highest observed exposure")
  expect_match(msg, "100")
})

test_that("an EC50 inside the exposure range gives no EC50 warning", {
  d   <- emax_data()
  res <- collect_warnings(fit_logistic_emax(d$y, d$x))
  expect_false(has_class(res$warnings, "bootlace_emax_ec50_warning"))
})

test_that("bootstrap_emax() keeps replicates that only fail identifiability checks", {
  set.seed(42)
  dose <- rep(c(0, 5, 10, 25, 50, 100), length.out = 200)
  y    <- rbinom(200, 1, plogis(qlogis(0.10) +
                                  (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)))
  fit  <- fit_logistic_emax(y, dose)
  res  <- collect_warnings(suppressMessages(
    bootstrap_emax(fit, n_boot = 120L, seed = 1L)
  ))
  # Replicate 25 of this run puts EC50 near 3000: kept, and reported.
  expect_equal(res$value$n_boot, 120L)
  expect_true(has_class(res$warnings, "bootlace_emax_ec50_warning"))
  msg <- conditionMessage(Filter(
    function(w) inherits(w, "bootlace_emax_ec50_warning"), res$warnings
  )[[1]])
  expect_match(msg, "of 120 bootstrap replicates")
})

# ---- Analytic gradient ---------------------------------------------------------

test_that(".emax_grad_fn() matches a numerical gradient of .emax_nll_fn()", {
  set.seed(11)
  n  <- 150
  x  <- rep(c(0, 2, 5, 10, 40), length.out = n)
  zl <- cbind(a = rnorm(n))
  ze <- cbind(b = rbinom(n, 1, 0.5))
  zc <- cbind(c = rnorm(n))
  y  <- rbinom(n, 1, 0.4)
  idx <- list(e0 = 1L, emax0 = 2L, log_ec50 = 3L, lin = 4L, emax_cov = 5L,
              ec50_cov = 6L)
  nll <- .emax_nll_fn(x, y, zl, ze, zc, idx)
  gr  <- .emax_grad_fn(x, y, zl, ze, zc, idx)
  for (th in list(c(-1, 2, log(5), 0.3, -0.7, 0.2),
                  c(0.5, -3, log(0.5), -1, 1.5, -0.4))) {
    expect_equal(gr(th), .num_grad(nll, th, h = 1e-6), tolerance = 1e-6)
  }
  # No covariates
  none <- matrix(0, n, 0L)
  idx0 <- list(e0 = 1L, emax0 = 2L, log_ec50 = 3L, lin = integer(0),
               emax_cov = integer(0), ec50_cov = integer(0))
  nll0 <- .emax_nll_fn(x, y, none, none, none, idx0)
  gr0  <- .emax_grad_fn(x, y, none, none, none, idx0)
  expect_equal(gr0(c(-1, 2, log(5))), .num_grad(nll0, c(-1, 2, log(5)), h = 1e-6),
               tolerance = 1e-6)
})
