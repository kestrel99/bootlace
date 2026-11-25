# Tests for the internal sample-size (EPV) helpers in R/sample_size.R.

make_bin <- function(n, p, seed = 1L) {
  set.seed(seed)
  data.frame(
    y  = rbinom(n, 1, p),
    x1 = rnorm(n),
    x2 = rnorm(n),
    g  = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
}

make_emax <- function(n = 200L, seed = 42L) {
  set.seed(seed)
  dose <- rep(c(0, 5, 10, 25, 50, 100), length.out = n)
  eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
  list(y = rbinom(n, 1, plogis(eta)), dose = dose)
}

# Minimal diagnostics list for testing .sample_size_problems() directly.
ss_stub <- function(epv, n_params = 4L) {
  list(
    type = "binary", limiting_events = epv * n_params, n_params = n_params,
    epv = epv, sparse_levels = character(0), n_x_levels = NA_integer_
  )
}

# ---- Counting ----------------------------------------------------------------

test_that(".sample_size() counts events and parameters for a binomial glm", {
  d   <- make_bin(200, 0.3)
  fit <- glm(y ~ x1 + x2 + g, data = d, family = binomial)
  ss  <- .sample_size(fit)
  expect_equal(ss$events, sum(d$y))
  expect_equal(ss$non_events, sum(1 - d$y))
  expect_equal(ss$limiting_events, min(sum(d$y), sum(1 - d$y)))
  expect_equal(ss$n_params, 4L)
  expect_equal(ss$epv, ss$limiting_events / 4)
})

test_that(".sample_size() weights a cbind() response by the number of trials", {
  agg <- data.frame(s = c(5, 12), f = c(15, 8), g = factor(c("a", "b")))
  fit <- glm(cbind(s, f) ~ g, data = agg, family = binomial)
  ss  <- .sample_size(fit)
  expect_equal(ss$events, 17)
  expect_equal(ss$non_events, 23)
  expect_equal(ss$n_params, 1L)
})

test_that(".sample_size() uses nevent for coxph and doesn't count strata", {
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + factor(sex) + strata(inst),
    data = survival::lung
  )
  ss <- .sample_size(fit)
  expect_equal(ss$events, fit$nevent)
  expect_true(is.na(ss$non_events))
  expect_equal(ss$n_params, 2L)
  expect_length(ss$sparse_levels, 0L)
})

test_that(".sample_size() counts Emax parameters excluding E0", {
  d   <- make_emax()
  fit <- suppressWarnings(fit_logistic_emax(d$y, d$dose))
  ss  <- .sample_size(fit)
  expect_equal(ss$n_params, 2L)
  expect_equal(ss$events, sum(d$y))
  expect_equal(ss$n_x_levels, 6L)
})

# ---- Checks ------------------------------------------------------------------

test_that("adequate data raises no problems and no warning", {
  d   <- make_bin(400, 0.4)
  fit <- glm(y ~ x1 + x2, data = d, family = binomial)
  ss  <- .sample_size(fit)
  expect_length(ss$problems, 0L)
  expect_no_warning(.warn_sample_size(ss, "f"))
})

test_that("low EPV and few events raise one classed warning", {
  d   <- make_bin(60, 0.2, seed = 3L)
  fit <- glm(y ~ x1 + x2 + g, data = d, family = binomial)
  ss  <- .sample_size(fit)
  expect_match(ss$problems, "is below 10", all = FALSE)
  expect_match(ss$problems, "fewer than 20", all = FALSE)
  expect_warning(
    .warn_sample_size(ss, "caller"),
    class = "bootlace_sample_size_warning"
  )
  expect_warning(.warn_sample_size(ss, "caller"), "^caller\\(\\)")
})

test_that("EPV of exactly 10 passes the fitting rule", {
  expect_length(.sample_size_problems(ss_stub(10)), 0L)
  expect_match(.sample_size_problems(ss_stub(9.9)), "is below 10")
})

test_that("selection tiers follow Heinze et al. (2018) Table 3", {
  expect_match(
    .sample_size_problems(ss_stub(10), selection = TRUE, alpha = 0.05),
    "not recommended"
  )
  mid <- .sample_size_problems(ss_stub(20), selection = TRUE, alpha = 0.05)
  expect_match(mid, "shrinkage")
  expect_match(mid, "0.157")
  mid_aic <- .sample_size_problems(ss_stub(20), selection = TRUE, alpha = 0.2)
  expect_false(grepl("0.157", mid_aic))
  expect_match(
    .sample_size_problems(ss_stub(25), selection = TRUE, alpha = 0.05),
    "at most 25"
  )
  expect_length(
    .sample_size_problems(ss_stub(26), selection = TRUE, alpha = 0.05),
    0L
  )
})

test_that("sparse factor levels are reported for glm", {
  d <- make_bin(300, 0.4)
  d$y[d$g == "c"] <- 0
  fit <- suppressWarnings(glm(y ~ x1 + g, data = d, family = binomial))
  ss  <- .sample_size(fit)
  expect_equal(ss$sparse_levels, "g = c (0 events)")
  expect_match(ss$problems, "Sparse factor levels", all = FALSE)
})

test_that("sparse factor levels are reported for coxph", {
  set.seed(2)
  d <- data.frame(
    time   = rexp(200),
    status = rbinom(200, 1, 0.7),
    g      = factor(sample(c("a", "b"), 200, replace = TRUE))
  )
  d$status[d$g == "b"] <- 0
  fit <- suppressWarnings(
    survival::coxph(survival::Surv(time, status) ~ g, data = d)
  )
  expect_equal(.sample_size(fit)$sparse_levels, "g = b (0 events)")
})

test_that("Emax sparse levels are rebuilt from factor covariates", {
  d <- make_emax()
  set.seed(1)
  grp <- factor(sample(c("A", "B", "C"), 200, replace = TRUE))
  y <- d$y
  y[grp == "C"] <- 0
  fit <- suppressWarnings(
    fit_logistic_emax(y, d$dose, emax_covs = data.frame(grp = grp))
  )
  expect_equal(.sample_size(fit)$sparse_levels, "grp = C (0 events)")
  fit_bare <- suppressWarnings(fit_logistic_emax(y, d$dose, emax_covs = grp))
  expect_equal(.sample_size(fit_bare)$sparse_levels, "emax_covs = C (0 events)")
})

test_that("binary 0/1 Emax covariate columns are checked", {
  d <- make_emax()
  set.seed(2)
  b <- rbinom(200, 1, 0.5)
  y <- d$y
  y[b == 1] <- 1
  fit <- suppressWarnings(
    fit_logistic_emax(y, d$dose, linear_covs = cbind(b = b))
  )
  expect_equal(.sample_size(fit)$sparse_levels, "b = 1 (0 non-events)")
})

test_that("Emax with fewer than 4 distinct exposure levels is flagged", {
  d    <- make_emax()
  dose <- rep(c(0, 10, 100), length.out = 200)
  fit  <- suppressWarnings(fit_logistic_emax(d$y, dose))
  ss   <- .sample_size(fit)
  expect_equal(ss$n_x_levels, 3L)
  expect_match(ss$problems, "exposure level", all = FALSE)
})

# ---- Print line --------------------------------------------------------------

test_that(".events_line() formats binary and survival counts", {
  expect_equal(
    .events_line(list(type = "binary", events = 54, non_events = 66,
                      n_params = 5L, epv = 10.8)),
    "54 events, 66 non-events (EPV 10.8, 5 params)"
  )
  expect_equal(
    .events_line(list(type = "survival", events = 165, n = 228,
                      n_params = 2L, epv = 82.5)),
    "165 of 228 (EPV 82.5, 2 params)"
  )
  expect_equal(
    .events_line(list(type = "binary", events = 10, non_events = 5,
                      n_params = 0L, epv = Inf)),
    "10 events, 5 non-events (EPV n/a, 0 params)"
  )
})

# ---- Entry points ------------------------------------------------------------

small_lung <- function() {
  d <- survival::lung[1:30, ]
  d$status <- d$status - 1L
  d
}

test_that("bootstrap_logistic() warns and records sample size for low EPV", {
  set.seed(5)
  d <- data.frame(y = rbinom(60, 1, 0.2), a = rnorm(60), b = rnorm(60),
                  c = rnorm(60))
  fit <- glm(y ~ a + b + c, data = d, family = binomial)
  expect_warning(
    bl <- suppressMessages(bootstrap_logistic(fit, n_boot = 150L, seed = 1L)),
    class = "bootlace_sample_size_warning"
  )
  expect_lt(bl$sample_size$epv, 10)
  expect_output(print(bl), "EPV")
})

test_that("bootstrap_cox() warns and records sample size for low EPV", {
  d   <- small_lung()
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + sex + ph.ecog, data = d
  )
  expect_warning(
    bc <- suppressMessages(bootstrap_cox(fit, n_boot = 150L, seed = 1L)),
    class = "bootlace_sample_size_warning"
  )
  expect_lt(bc$sample_size$epv, 10)
  expect_output(print(bc), "EPV")
})

test_that("fit_logistic_emax() warns and records sample size", {
  d    <- make_emax()
  dose <- rep(c(0, 10, 100), length.out = 200)
  expect_warning(
    fit <- fit_logistic_emax(d$y, dose),
    class = "bootlace_sample_size_warning"
  )
  expect_equal(fit$sample_size$n_x_levels, 3L)
  expect_output(print(fit), "EPV")
})

test_that("bootstrap_emax() doesn't count sample-size warnings as failed refits", {
  d   <- make_emax()
  fit <- fit_logistic_emax(d$y, d$dose)
  # Compare with the same bootstrap run without the extra warning. The EC50
  # range warning is expected here and is tested in test-fit_logistic_emax.R.
  plain <- suppressMessages(muffle_class(
    bootstrap_emax(fit, n_boot = 120L, seed = 1L), "bootlace_emax_ec50_warning"
  ))
  real_fit <- fit_logistic_emax
  local_mocked_bindings(fit_logistic_emax = function(...) {
    .warn_sample_size(list(problems = "test problem"), "fit_logistic_emax")
    real_fit(...)
  })
  be <- suppressMessages(muffle_class(
    bootstrap_emax(fit, n_boot = 120L, seed = 1L), "bootlace_emax_ec50_warning"
  ))
  expect_equal(be$n_boot, plain$n_boot)
  expect_gt(be$n_boot, 100L)
  expect_output(print(be), "EPV")
})

test_that("stepwise_reduce() applies the selection tiers to the full model", {
  d   <- make_bin(200, 0.3, seed = 4L)   # 68 events, 4 parameters: EPV 17
  fit <- glm(y ~ x1 + x2 + g, data = d, family = binomial)
  expect_warning(
    sr <- stepwise_reduce(fit, alpha = 0.05, verbose = FALSE),
    "0.157",
    class = "bootlace_sample_size_warning"
  )
  expect_equal(sr$sample_size$epv_threshold, 25)
  expect_output(print(sr), "EPV")
})

test_that("stepwise_reduce() says selection is not recommended at EPV <= 10", {
  d   <- make_bin(80, 0.3, seed = 4L)
  fit <- glm(y ~ x1 + x2 + g, data = d, family = binomial)
  expect_warning(
    stepwise_reduce(fit, verbose = FALSE),
    "not recommended",
    class = "bootlace_sample_size_warning"
  )
})

test_that("vpc_cox() warns and records sample size for a small Cox fit", {
  d   <- small_lung()
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + sex + ph.ecog, data = d
  )
  expect_warning(
    vpc <- vpc_cox(fit, d, n_sim = 5, n_grid = 10, seed = 1),
    class = "bootlace_sample_size_warning"
  )
  expect_lt(vpc$sample_size$epv, 10)
})
