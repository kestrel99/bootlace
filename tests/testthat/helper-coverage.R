# Simulation helpers for the opt-in coverage study in test-coverage.R
# (SCIENTIFIC_REVIEW.md "Verification priorities"). Each scenario simulates
# data with known parameters, runs the package bootstrap, and returns the
# proportion of simulations in which each 95% percentile interval covers the
# truth.

COVERAGE_N_SIM  <- 100L
COVERAGE_N_BOOT <- 150L

# Run `expr` silently: muffle messages and every warning.
.quietly <- function(expr) {
  suppressMessages(withCallingHandlers(
    expr, warning = function(w) invokeRestart("muffleWarning")
  ))
}

# TRUE/FALSE per column: does the percentile interval of `boot` cover `truth`?
.covers <- function(boot, truth, level = 0.95) {
  a <- (1 - level) / 2
  q <- apply(boot, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE)
  q[1, ] <= truth & truth <= q[2, ]
}

# Proportion of simulations covering the truth. `one_sim(s)` returns either a
# bootstrap matrix (percentile intervals are checked against `truth`) or a
# named logical vector of hits.
.coverage <- function(one_sim, truth, n_sim = COVERAGE_N_SIM) {
  hits <- lapply(seq_len(n_sim), function(s) {
    set.seed(1000L + s)
    res <- one_sim(s)
    if (is.matrix(res)) .covers(res, truth) else res
  })
  rowMeans(do.call(cbind, hits))
}

coverage_logistic <- function(n = 300L) {
  truth <- c("(Intercept)" = -1, x1 = 0.8, x2 = -0.5)
  .coverage(function(s) {
    d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rbinom(n, 1, 0.5))
    d$y <- stats::rbinom(n, 1, stats::plogis(truth[[1]] + truth[[2]] * d$x1 +
                                               truth[[3]] * d$x2))
    fit <- stats::glm(y ~ x1 + x2, data = d, family = stats::binomial)
    .quietly(bootstrap_logistic(fit, n_boot = COVERAGE_N_BOOT, seed = s))$boot_coefs
  }, truth)
}

# Exponential event times with independent uniform censoring (about half of
# subjects have an event).
.cox_data <- function(n, truth) {
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rbinom(n, 1, 0.5))
  t_event <- stats::rexp(n, 0.1 * exp(truth[["x1"]] * d$x1 + truth[["x2"]] * d$x2))
  t_cens  <- stats::runif(n, 0, 25)
  d$time   <- pmin(t_event, t_cens)
  d$status <- as.integer(t_event <= t_cens)
  d
}

coverage_cox <- function(n = 300L) {
  truth <- c(x1 = 0.5, x2 = -0.7)
  .coverage(function(s) {
    d   <- .cox_data(n, truth)
    fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = d)
    .quietly(bootstrap_cox(fit, n_boot = COVERAGE_N_BOOT, seed = s))$boot_coefs
  }, truth)
}

# Patient-level resampling: each subject's follow-up is split into several
# (start, stop] rows; bootstrap_cox() must resample whole subjects.
coverage_cox_id <- function(n = 300L) {
  truth <- c(x1 = 0.5, x2 = -0.7)
  .coverage(function(s) {
    d <- .cox_data(n, truth)
    d$id <- seq_len(n)
    long <- survival::survSplit(data = d, cut = c(2, 5, 10), end = "time",
                                event = "status", start = "tstart")
    fit <- survival::coxph(survival::Surv(tstart, time, status) ~ x1 + x2,
                           data = long, id = id)
    .quietly(bootstrap_cox(fit, n_boot = COVERAGE_N_BOOT, seed = s))$boot_coefs
  }, truth)
}

# Emax: coverage of the intervals the package reports, BCa bootstrap
# (summary.boot_emax() default) and profile likelihood (confint() default).
.coverage_emax <- function(ec50, n = 300L) {
  truth <- c(e0 = stats::qlogis(0.1), emax = 3, log_ec50 = log(ec50))
  truth_rep <- c(e0 = truth[["e0"]], emax = 3, ec50 = ec50)
  doses <- c(0, 5, 10, 25, 50, 100)
  .coverage(function(s) {
    x <- rep(doses, length.out = n)
    y <- stats::rbinom(n, 1, stats::plogis(truth[["e0"]] +
                                             truth[["emax"]] * x / (ec50 + x)))
    fit <- .quietly(fit_logistic_emax(y, x))
    be  <- .quietly(bootstrap_emax(fit, n_boot = COVERAGE_N_BOOT, seed = s))
    bca <- vapply(1:3, function(j) {
      ci <- .bca_ci(be$boot_theta[, j], fit$theta_hat[[j]],
                    be$jack_theta[, j], 0.95)
      ci[1] <= truth[[j]] && truth[[j]] <= ci[2]
    }, logical(1L))
    pr  <- .quietly(confint(fit, parm = names(truth_rep), method = "profile"))
    prof <- pr[, 1] <= truth_rep & truth_rep <= pr[, 2]
    c(stats::setNames(bca, paste("BCa", names(truth))),
      stats::setNames(prof, paste("profile", names(truth_rep))))
  }, truth)
}

coverage_emax_identified <- function() .coverage_emax(ec50 = 10)
coverage_emax_weak       <- function() .coverage_emax(ec50 = 60)
