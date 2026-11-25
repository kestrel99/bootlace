# vpc_cox() rejects unsupported fits and inputs

    Code
      vpc_cox(lm(time ~ age, lung), lung)
    Condition
      Error:
      ! `fit` must be a coxph() fit
    Code
      vpc_cox(strat, lung)
    Condition
      Error:
      ! models with strata() aren't supported
    Code
      vpc_cox(counting, lung)
    Condition
      Error:
      ! only right-censored data are supported, not (start, stop] data
    Code
      vpc_cox(fit, lung[-1, ])
    Condition
      Error:
      ! `data` doesn't match the data `fit` was fitted to
    Code
      vpc_cox(fit, lung, by = "nope")
    Condition
      Error:
      ! `data` has no column nope
    Code
      vpc_cox(fit, lung, by = "many", overlay = TRUE)
    Condition
      Error:
      ! `overlay` supports up to 8 groups; use panels instead

