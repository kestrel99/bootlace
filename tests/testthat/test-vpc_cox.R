lung_df <- function() {
  lung <- survival::lung
  lung$status <- lung$status - 1L
  lung
}

lung_fit <- function(data = lung_df()) {
  survival::coxph(survival::Surv(time, status) ~ age + sex, data = data)
}

test_that("vpc_cox() returns simulated and observed curves and a plot", {
  vpc <- vpc_cox(lung_fit(), lung_df(), n_sim = 20, n_grid = 25, seed = 1)

  expect_named(vpc, c("sim", "obs", "zph", "sample_size", "plot"))
  expect_named(vpc$sim, c("group", "time", "lo", "median", "hi"))
  expect_equal(nrow(vpc$sim), 25)
  expect_equal(levels(vpc$sim$group), "All (n = 228)")
  expect_gte(min(vpc$sim$median - vpc$sim$lo), 0)
  expect_gte(min(vpc$sim$hi - vpc$sim$median), 0)
  expect_s3_class(vpc$zph, "cox.zph")
  expect_s3_class(vpc$plot, "ggplot")
})

test_that("vpc_cox() simulated median tracks the observed curve", {
  vpc <- vpc_cox(lung_fit(), lung_df(), n_sim = 100, n_grid = 50, seed = 1)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = lung_df())
  obs_at_grid <- summary(km, times = vpc$sim$time, extend = TRUE)$surv

  expect_lt(max(abs(vpc$sim$median - obs_at_grid)), 0.1)
})

test_that("vpc_cox() is reproducible with a seed", {
  run <- function() {
    vpc_cox(lung_fit(), lung_df(), n_sim = 10, n_grid = 10, seed = 3)$sim
  }
  expect_equal(run(), run())
})

test_that("vpc_cox() stratifies by a column, in panels or overlaid", {
  lung <- lung_df()
  for (overlay in c(FALSE, TRUE)) {
    vpc <- vpc_cox(
      lung_fit(lung),
      lung,
      by = "sex",
      n_sim = 10,
      n_grid = 10,
      seed = 1,
      overlay = overlay
    )
    expect_equal(levels(vpc$sim$group), c("1 (n = 138)", "2 (n = 90)"))
    expect_equal(nrow(vpc$sim), 20)
    expect_match(vpc$plot$labels$subtitle, "sex p = ")
  }
})

test_that("vpc_cox() drops rows the fit dropped for missing values", {
  lung <- lung_df()
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + ph.ecog,
    data = lung
  )
  vpc <- vpc_cox(fit, lung, n_sim = 5, n_grid = 10, seed = 1)
  expect_equal(levels(vpc$sim$group), "All (n = 227)")
})

test_that("vpc_cox() ends the observed curve at t_max", {
  vpc <- vpc_cox(
    lung_fit(),
    lung_df(),
    n_sim = 5,
    n_grid = 10,
    t_max = 500,
    seed = 1
  )
  expect_equal(max(vpc$obs$time), 500)
  expect_equal(max(vpc$sim$time), 500)
})

test_that("vpc_cox() supports bootstrap and reverse-KM censoring", {
  vpc <- vpc_cox(
    lung_fit(),
    lung_df(),
    n_sim = 10,
    n_grid = 10,
    seed = 1,
    bootstrap = TRUE,
    censoring = "reverseKM"
  )
  expect_equal(nrow(vpc$sim), 10)
  expect_match(vpc$plot$labels$caption, "bootstrap")
})

test_that(".reverse_km_sampler() draws censoring times after the event", {
  lung <- lung_df()
  draw <- .reverse_km_sampler(lung$time, lung$status, max(lung$time))
  set.seed(1)
  c_time <- draw()

  expect_length(c_time, sum(lung$status == 1))
  expect_gt(min(c_time - lung$time[lung$status == 1]), 0)
})

test_that("vpc_cox() rejects unsupported fits and inputs", {
  lung <- lung_df()
  fit <- lung_fit(lung)
  strat <- survival::coxph(
    survival::Surv(time, status) ~ age + survival::strata(sex),
    data = lung
  )
  lung$start <- 0
  counting <- survival::coxph(
    survival::Surv(start, time, status) ~ age,
    data = lung
  )
  lung$many <- rep(1:9, length.out = nrow(lung))

  expect_snapshot(error = TRUE, {
    vpc_cox(lm(time ~ age, lung), lung)
    vpc_cox(strat, lung)
    vpc_cox(counting, lung)
    vpc_cox(fit, lung[-1, ])
    vpc_cox(fit, lung, by = "nope")
    vpc_cox(fit, lung, by = "many", overlay = TRUE)
  })
})
