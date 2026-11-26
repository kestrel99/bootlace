## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* Local: Windows 11, R 4.6.1

## Notes for reviewers

* The examples for `bootstrap_emax()` and `plot_emax_curve()` take about
  12 seconds each, because each bootstrap replicate refits a non-linear
  model, so they are wrapped in `\donttest{}`. The model fit itself runs
  unconditionally.
* Method references are given in the Description field (Efron and Tibshirani
  1993; Peduzzi et al. 1996) and in the help pages.
