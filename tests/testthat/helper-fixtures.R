# Helper fixtures used across the testthat suite.
# Sourced automatically by testthat because the file name starts with "helper".

# ---- Small well-conditioned dataset (n = 120) --------------------------------
make_df <- function(n = 120, seed = 1L) {
  set.seed(seed)
  data.frame(
    age = rnorm(n, 50, 10),
    bmi = rnorm(n, 27, 5),
    sex = factor(sample(c("F", "M"), n, replace = TRUE)),
    trt = factor(
      sample(c("Placebo", "Low", "High"), n, replace = TRUE),
      levels = c("Placebo", "Low", "High")
    ),
    y = rbinom(n, 1, 0.45)
  )
}

DF <- make_df()
FIT <- glm(y ~ age + bmi + sex + trt, data = DF, family = binomial)

# Pre-built boot_logistic objects (low n_boot for speed)
BL <- bootstrap_logistic(FIT, n_boot = 150L, seed = 7L)

# Mute bootlace_sample_size_warning in tests that aren't about sample size,
# so other expectations about warnings still see only their own.
quiet_sample_size <- function(expr) {
  withCallingHandlers(
    expr,
    bootlace_sample_size_warning = function(w) invokeRestart("muffleWarning")
  )
}

# Mute warnings of the given condition class(es) only.
muffle_class <- function(expr, class) {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, class)) invokeRestart("muffleWarning")
  })
}
