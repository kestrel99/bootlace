# Coverage of 95% bootstrap percentile intervals in simulations with known
# parameters (SCIENTIFIC_REVIEW.md "Verification priorities").
#
# Opt-in: these take about 45 minutes. Run with
#   Sys.setenv(BOOTLACE_COVERAGE = "true"); devtools::test(filter = "coverage")
#
# With COVERAGE_N_SIM = 100 simulations the Monte Carlo SE of a coverage near
# 0.95 is about 0.022, so each interval must cover in at least 88% of
# simulations (about 3 SE below nominal). Simulation code: helper-coverage.R.
#
# Emax is checked with the intervals the package reports: BCa bootstrap
# (summary.boot_emax()) and profile likelihood (confint()). Plain percentile
# intervals under-covered Emax (0.86 well identified; 0.80-0.83 weakly
# identified), which is why those are the defaults. Results on 2026-11-25:
#   logistic                    0.97 0.93 0.95
#   Cox                         0.98 0.90
#   Cox, patient-level          0.98 0.90
#   Emax, EC50 = 10, BCa        0.90 0.93 0.95   profile 0.92 0.95 0.94
#   Emax, EC50 = 60, BCa        0.89 0.90 0.93   profile 0.91 0.99 0.95
# The weakly identified design also meets the 0.88 bar, so it shares it.

skip_if_not(
  identical(Sys.getenv("BOOTLACE_COVERAGE"), "true"),
  "Set BOOTLACE_COVERAGE=true to run the coverage simulations."
)

min_coverage <- 0.88

expect_coverage <- function(cov) {
  # Opt-in study: show the coverage values, not just pass/fail.
  cat(sprintf("\n  coverage: %s\n",
              paste(sprintf("%s %.2f", names(cov), cov), collapse = ", ")))
  for (nm in names(cov))
    expect_gte(cov[[nm]], min_coverage,
               label = sprintf("coverage of %s (%.2f)", nm, cov[[nm]]))
}

test_that("bootstrap_logistic() percentile intervals cover", {
  expect_coverage(coverage_logistic())
})

test_that("bootstrap_cox() percentile intervals cover", {
  expect_coverage(coverage_cox())
})

test_that("bootstrap_cox() with patient-level resampling covers", {
  expect_coverage(coverage_cox_id())
})

test_that("Emax BCa and profile intervals cover in a well-identified design", {
  expect_coverage(coverage_emax_identified())
})

test_that("Emax BCa and profile intervals cover in a weakly identified design", {
  expect_coverage(coverage_emax_weak())
})
