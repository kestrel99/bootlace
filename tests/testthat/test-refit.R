# Regression tests for SCIENTIFIC_REVIEW.md finding 1: bootstrap refits and
# stepwise LRTs must keep the original analysis rows, weights, offsets and
# (for Cox) strata and clusters.

library(survival)

make_pop <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    age   = rnorm(n, 50, 10),
    noise = rnorm(n),
    grp   = factor(sample(c("a", "b"), n, replace = TRUE)),
    w     = sample(1:3, n, replace = TRUE),
    expo  = log(runif(n, 1, 3))
  )
  d$y <- rbinom(n, 1, plogis(-3 + 0.05 * d$age))
  d$noise[sample(n, 100)] <- NA        # a third of rows missing a predictor
  d
}

# ---- GLM design and refit ------------------------------------------------------

test_that(".refit_glm() on all rows reproduces a weighted, offset, subset fit", {
  d   <- make_pop()
  fit <- glm(y ~ age + grp + noise, data = d, family = binomial,
             weights = w, offset = expo, subset = age > 35)
  des <- .glm_design(fit)
  b   <- .refit_glm(des, seq_len(des$n))
  expect_equal(b, coef(fit), tolerance = 1e-8)
})

test_that(".analysis_data() returns only the rows the fit used", {
  d   <- make_pop()
  fit <- glm(y ~ age + noise, data = d, family = binomial, subset = age > 35)
  ad  <- .analysis_data(fit)
  expect_equal(nrow(ad), nrow(model.frame(fit)))
  expect_false(anyNA(ad$noise))
  expect_true(all(ad$age > 35))
})

test_that("bootstrap_logistic() resamples only the analysis rows", {
  d   <- make_pop()
  fit <- glm(y ~ age + noise, data = d, family = binomial)
  bl  <- quiet_sample_size(bootstrap_logistic(fit, n_boot = 100L, seed = 1L))
  expect_equal(nrow(bl$data), nrow(model.frame(fit)))
})

test_that("bootstrap_logistic() keeps case weights in every refit", {
  d     <- make_pop()
  d$w2  <- ifelse(d$y == 1, 5, 1)      # strong weights shift the intercept
  fit_w <- glm(y ~ age, data = d, family = binomial, weights = w2)
  bl    <- quiet_sample_size(suppressWarnings(
    bootstrap_logistic(fit_w, n_boot = 200L, seed = 1L)
  ))
  # The bootstrap distribution must be centred on the weighted estimate,
  # not on the unweighted one.
  fit_u <- glm(y ~ age, data = d, family = binomial)
  centre <- median(bl$boot_coefs[, "(Intercept)"])
  expect_lt(abs(centre - coef(fit_w)[[1]]), abs(centre - coef(fit_u)[[1]]))
})

# ---- Cox design and refit ------------------------------------------------------

test_that(".refit_cox() on all rows reproduces a weighted, stratified, offset fit", {
  lung <- survival::lung
  lung$w <- runif(nrow(lung), 0.5, 2)
  fit <- coxph(Surv(time, status) ~ age + ph.ecog + offset(0.01 * meal.cal) +
                 strata(sex), data = lung, weights = w)
  des <- .cox_design(fit)
  b   <- .refit_cox(des, seq_len(des$n))
  expect_equal(b, coef(fit), tolerance = 1e-6)
})

test_that(".cox_design() picks up the subject identifier", {
  fit_id <- coxph(Surv(tstart, tstop, status) ~ age + sex,
                  data = survival::cgd, id = id)
  expect_false(is.null(.cox_design(fit_id)$cluster))
  fit_cl <- coxph(Surv(time, status) ~ age + cluster(inst), data = survival::lung)
  expect_false(is.null(.cox_design(fit_cl)$cluster))
})

test_that(".resample_rows() keeps every row of a resampled cluster together", {
  des <- list(n = 6L, cluster = c("a", "b", "a", "c", "b", "c"))
  set.seed(1)
  rows <- .resample_rows(des)
  # 3 clusters drawn, each contributing both of its rows, consecutively
  expect_length(rows, 6L)
  pairs <- matrix(des$cluster[rows], ncol = 2, byrow = TRUE)
  expect_true(all(pairs[, 1] == pairs[, 2]))
  # plain row resampling when there are no clusters
  expect_length(.resample_rows(list(n = 6L, cluster = NULL)), 6L)
})

test_that("bootstrap_cox() rejects counting-process data without an id", {
  fit <- coxph(Surv(tstart, tstop, status) ~ age + sex, data = survival::cgd)
  expect_error(bootstrap_cox(fit, n_boot = 10L), "id")
})

# ---- Stepwise ------------------------------------------------------------------

test_that("stepwise_reduce() LRTs use the full model's rows", {
  d   <- make_pop()
  fit <- glm(y ~ age + noise, data = d, family = binomial)
  sr  <- quiet_sample_size(stepwise_reduce(fit, alpha = 0.05, verbose = FALSE))
  ref <- drop1(fit, test = "LRT")       # drop1 keeps the model frame's rows
  expect_equal(sr$step_table$d_neg2LL[1], ref["noise", "LRT"], tolerance = 1e-6)
})

test_that("stepwise_reduce() keeps case weights in candidate refits", {
  d   <- make_pop()
  d   <- d[!is.na(d$noise), ]
  # Non-integer weights: logLik() of a binomial GLM is then not a likelihood,
  # so the LRT has to come from the deviance, as drop1() does.
  d$w <- runif(nrow(d), 0.5, 2)
  fit <- suppressWarnings(
    glm(y ~ age + noise, data = d, family = binomial, weights = w)
  )
  sr  <- quiet_sample_size(suppressWarnings(
    stepwise_reduce(fit, alpha = 0.05, verbose = FALSE)
  ))
  ref <- suppressWarnings(drop1(fit, test = "LRT"))
  expect_equal(sr$step_table$d_neg2LL[1], ref["noise", "LRT"], tolerance = 1e-6)
})
