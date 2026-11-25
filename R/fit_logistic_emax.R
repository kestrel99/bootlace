#' Fit a logistic Emax model
#'
#' Fits a logistic Emax model to binary outcome data via maximum likelihood,
#' with optional additive linear covariates and covariate effects on the Emax
#' and EC50 parameters.
#'
#' The model for observation \eqn{i} is
#' \deqn{
#'   \eta_i = E_0
#'     + \bigl(E_{max} + \mathbf{z}_{emax,i}^\top\boldsymbol{\gamma}\bigr)
#'       \frac{x_i}{e^{\ell_i} + x_i}
#'     + \mathbf{z}_{lin,i}^\top\boldsymbol{\beta},
#' }
#' where \eqn{\ell_i = \log(EC_{50}) + \mathbf{z}_{ec50,i}^\top\boldsymbol{\delta}}
#' keeps \eqn{EC_{50,i} > 0}, and
#' \eqn{P(Y_i=1) = \text{logistic}(\eta_i)}.
#' When all covariate arguments are `NULL` the model reduces to the standard
#' three-parameter logistic Emax formula.
#'
#' @section Covariate matrices:
#' `linear_covs`, `emax_covs`, and `ec50_covs` accept numeric matrices,
#' factors, character vectors, or data frames (possibly mixing numeric and
#' categorical columns).  Categorical variables are automatically expanded via
#' [model.matrix()] using treatment contrasts: for a factor with \eqn{k}
#' levels the reference level is absorbed by the base parameter (\eqn{E_0},
#' \eqn{E_{max}}, or \eqn{\log EC_{50}}) and \eqn{k-1} contrast columns are
#' added.  All inputs must have `length(y)` rows.  Column names label the
#' coefficients; unnamed columns receive labels `"V1"`, `"V2"`, etc.
#' When predicting with categorical covariates, pass new data in the same form
#' (factor or data frame); [predict.logistic_emax()] re-applies the training
#' contrasts and reference levels automatically.
#'
#' @section Exposure levels:
#' `x` must be finite and non-negative, and `y` must contain both outcome
#' classes. E0, Emax and EC50 need at least 3 distinct exposure values
#' (counting zero), so fewer is an error. With exactly 3 the curve is exactly
#' determined and EC50 is usually poorly estimated, so in addition to the
#' **Sample size checks** a warning is given when `x` has fewer than 4
#' distinct values. For factor covariates, the sparse-level check uses the
#' original levels. The EPV count is all parameters except E0. EPV rules of
#' thumb were developed for models with a linear predictor. Emax and EC50
#' enter this model nonlinearly, so treat EPV >= 10 as a minimum.
#'
#' @section Convergence and standard errors:
#' When `start` is `NULL`, the optimiser is run from three EC50 starting
#' values (the 25th, 50th and 75th percentiles of the non-zero exposures) and
#' the best optimum is kept. The result is then finished with Newton steps.
#' `converged` is `TRUE` only when [optim()] reports success, the Hessian at
#' the optimum is positive definite and well conditioned, and the remaining
#' possible gain in log-likelihood (the Newton decrement) is below 1e-4.
#' Otherwise a warning explains which check failed. When the Hessian check
#' fails, some parameters are not identified (for example EC50 far above the
#' highest exposure, where Emax and EC50 trade off), and standard errors are
#' reported as `NA` instead of as misleadingly precise values. In weakly
#' identified designs, consider a profile likelihood or [bootstrap_emax()]
#' for EC50 rather than the Wald standard error.
#'
#' @section Unrealistic EC50:
#' A detailed warning of class `bootlace_emax_ec50_warning` is given when the
#' EC50 estimate (for the reference covariate values) lies outside the range
#' the observed exposures can support: above the highest exposure, where the
#' plateau is never observed and Emax and EC50 trade off against each other,
#' or below a tenth of the lowest non-zero exposure, where the curve is
#' effectively a step. The warning explains the consequence and lists
#' options. The limits are this package's own rule of thumb.
#'
#' The identifiability and stopping warnings described under **Convergence
#' and standard errors** have class `bootlace_emax_identification_warning`.
#'
#' @param y Numeric binary outcome vector (0 or 1).
#' @param x Numeric predictor vector (e.g. dose, concentration, time), same
#'   length as `y`.
#' @param linear_covs Optional covariates that enter the linear predictor
#'   additively (\eqn{\mathbf{z}_{lin}}).  Accepts a numeric matrix, a factor,
#'   a character vector, or a data frame; see **Covariate matrices**. `NULL`
#'   omits.
#' @param emax_covs Optional covariates that modify the Emax parameter
#'   additively (\eqn{\mathbf{z}_{emax}}).  Same formats as `linear_covs`.
#'   `NULL` omits.
#' @param ec50_covs Optional covariates that modify \eqn{\log(EC_{50})}
#'   additively (\eqn{\mathbf{z}_{ec50}}); a coefficient of \eqn{\delta}
#'   multiplies EC50 by \eqn{e^\delta} per unit change.  Same formats as
#'   `linear_covs`.  `NULL` omits.
#' @param start Starting values, which must be finite with a positive EC50.
#'   One of:
#'   \itemize{
#'     \item `NULL` (default): auto-derived from a simple logistic regression,
#'       with three EC50 starting values (see **Convergence and standard
#'       errors**); covariate effects start at 0.
#'     \item A length-3 vector `c(e0, emax, ec50)` on the natural scale:
#'       base parameters; covariate effects start at 0.
#'     \item A full-length vector of length
#'       `3 + ncol(linear_covs) + ncol(emax_covs) + ncol(ec50_covs)` ordered
#'       `(e0, emax, ec50, beta..., gamma..., delta...)`, all on the natural
#'       scale (ec50 is log-transformed internally).
#'   }
#' @param control A list of control parameters passed to [optim()]. Default
#'   `list(maxit = 1000)`.
#' @param hessian Logical; if `TRUE` (default) the Hessian is computed and used
#'   to derive the variance-covariance matrix and standard errors.
#' @param ... Additional arguments passed to [optim()].
#'
#' @return An S3 object of class `"logistic_emax"`, a named list with:
#'   \describe{
#'     \item{`coefficients`}{Named numeric vector on the reporting scale:
#'       `c(e0, emax, ec50, linear.*, emax.*, log_ec50.*)`.  `ec50` is
#'       back-transformed; EC50 covariate effects (`log_ec50.*`) remain on the
#'       log scale.}
#'     \item{`std_err`}{Standard errors matching `coefficients` (delta method
#'       applied to `ec50`).}
#'     \item{`vcov_theta`}{Variance-covariance matrix in the optimisation
#'       parameterisation `(e0, emax0, log_ec50_0, beta..., gamma...,
#'       delta...)`, or `NULL` when `hessian = FALSE` or the Hessian is not
#'       positive definite.}
#'     \item{`theta_hat`}{Raw optimised parameter vector (same parameterisation
#'       as `vcov_theta`).}
#'     \item{`fitted.values`}{Fitted probabilities for the observed data.}
#'     \item{`linear.predictors`}{Fitted values on the logit scale.}
#'     \item{`residuals`}{Raw residuals `y - fitted.values`.}
#'     \item{`loglik`}{Log-likelihood at the MLE.}
#'     \item{`AIC`}{Akaike information criterion.}
#'     \item{`converged`}{`TRUE` when the fit passed every convergence check;
#'       see **Convergence and standard errors**.}
#'     \item{`n_par`}{Total number of parameters.}
#'     \item{`x`, `y`}{Complete-case predictor and outcome vectors.}
#'     \item{`Z_lin`, `Z_emax`, `Z_ec50`}{Covariate matrices used in fitting
#'       (0-column matrices when the corresponding argument was `NULL`).}
#'     \item{`idx`}{Named list of integer indices into `theta_hat` for each
#'       parameter group, used internally by [predict.logistic_emax()].}
#'     \item{`mm_info_lin`, `mm_info_emax`, `mm_info_ec50`}{List (or `NULL`)
#'       of model-matrix metadata (terms, xlevels, contrasts) present when
#'       categorical covariates were supplied; consumed by
#'       [predict.logistic_emax()] to re-expand new data consistently.}
#'     \item{`call`}{The matched function call.}
#'     \item{`sample_size`}{Sample-size diagnostics; see
#'       **Sample size checks**.}
#'   }
#'
#' @seealso [predict.logistic_emax()], [print.logistic_emax()],
#'   [roc_logistic_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#'
#' # Basic model
#' fit <- fit_logistic_emax(y, dose)
#' print(fit)
#'
#' # Linear covariate on the link scale
#' wt  <- rnorm(n, 70, 10)
#' fit2 <- fit_logistic_emax(y, dose, linear_covs = cbind(weight = wt))
#' print(fit2)
#'
#' # Covariate modifying Emax (binary group as numeric)
#' grp <- rbinom(n, 1, 0.5)
#' fit3 <- fit_logistic_emax(y, dose, emax_covs = cbind(group = grp))
#' print(fit3)
#'
#' # Categorical covariate: factor passed directly (no pre-expansion needed)
#' trt <- factor(sample(c("placebo", "low", "high"), n, replace = TRUE))
#' fit4 <- fit_logistic_emax(y, dose, linear_covs = data.frame(treatment = trt))
#' print(fit4)
#'
#' @template sample-size
#' @export
fit_logistic_emax <- function(
  y,
  x,
  linear_covs = NULL,
  emax_covs   = NULL,
  ec50_covs   = NULL,
  start       = NULL,
  control     = list(maxit = 1000),
  hessian     = TRUE,
  ...
) {
  if (length(y) != length(x)) stop("y and x must have the same length.")
  if (!all(y %in% c(0, 1)))   stop("y must be binary (0/1).")

  # ---- Covariate matrices (original length, before subsetting) ---------------
  to_mat <- function(z, n_orig, label) {
    if (is.null(z)) return(list(mat = matrix(0, nrow = n_orig, ncol = 0), mm_info = NULL))
    needs_expansion <- is.factor(z) || is.character(z) ||
      (is.data.frame(z) &&
         any(vapply(z, function(v) is.factor(v) || is.character(v), logical(1L))))
    if (needs_expansion) {
      df  <- if (is.data.frame(z)) z else
        data.frame(.var = if (is.character(z)) as.factor(z) else z)
      df[]  <- lapply(df, function(v) if (is.character(v)) as.factor(v) else v)
      mf  <- stats::model.frame(~ ., data = df, na.action = stats::na.pass)
      mm  <- stats::model.matrix(attr(mf, "terms"), data = mf)
      mm  <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      xlevels <- lapply(mf, function(v) if (is.factor(v)) levels(v) else NULL)
      xlevels <- xlevels[!vapply(xlevels, is.null, logical(1L))]
      mm_info <- list(
        terms     = attr(mf, "terms"),
        xlevels   = xlevels,
        contrasts = attr(mm, "contrasts")
      )
      m <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm), dimnames = dimnames(mm))
    } else {
      m       <- as.matrix(z)
      if (!is.numeric(m))
        stop(label, " must be numeric, a factor, or a data frame.")
      mm_info <- NULL
    }
    if (nrow(m) != n_orig) stop(label, " must have ", n_orig, " rows.")
    if (ncol(m) > 0L && is.null(colnames(m)))
      colnames(m) <- paste0("V", seq_len(ncol(m)))
    list(mat = m, mm_info = mm_info)
  }
  n_orig       <- length(y)
  res_lin      <- to_mat(linear_covs, n_orig, "linear_covs")
  res_emax     <- to_mat(emax_covs,   n_orig, "emax_covs")
  res_ec50     <- to_mat(ec50_covs,   n_orig, "ec50_covs")
  Z_lin_raw    <- res_lin$mat;   mm_info_lin  <- res_lin$mm_info
  Z_emax_raw   <- res_emax$mat;  mm_info_emax <- res_emax$mm_info
  Z_ec50_raw   <- res_ec50$mat;  mm_info_ec50 <- res_ec50$mm_info

  # ---- Complete-case subsetting ----------------------------------------------
  ok <- stats::complete.cases(y, x)
  if (ncol(Z_lin_raw)  > 0) ok <- ok & stats::complete.cases(Z_lin_raw)
  if (ncol(Z_emax_raw) > 0) ok <- ok & stats::complete.cases(Z_emax_raw)
  if (ncol(Z_ec50_raw) > 0) ok <- ok & stats::complete.cases(Z_ec50_raw)

  y      <- y[ok]; x <- x[ok]; n <- length(y)
  if (any(!is.finite(x))) stop("`x` must be finite.")
  if (any(x < 0))         stop("`x` must be non-negative.")
  if (length(unique(y)) < 2L)
    stop("`y` must contain both outcome classes (0 and 1).")
  n_x <- length(unique(x))
  if (n_x < 3L)
    stop(sprintf(paste0(
      "`x` has %d distinct value(s); at least 3 (counting zero) are needed ",
      "to estimate E0, Emax and EC50."
    ), n_x))
  Z_lin  <- Z_lin_raw[ok,  , drop = FALSE]
  Z_emax <- Z_emax_raw[ok, , drop = FALSE]
  Z_ec50 <- Z_ec50_raw[ok, , drop = FALSE]

  p_lin  <- ncol(Z_lin)
  p_emax <- ncol(Z_emax)
  p_ec50 <- ncol(Z_ec50)
  n_par  <- 3L + p_lin + p_emax + p_ec50

  # ---- Parameter indices (1-based into theta) --------------------------------
  idx_e0       <- 1L
  idx_emax0    <- 2L
  idx_log_ec50 <- 3L
  idx_lin      <- if (p_lin  > 0) seq.int(4L,             3L + p_lin)                    else integer(0)
  idx_emax_cov <- if (p_emax > 0) seq.int(4L + p_lin,     3L + p_lin + p_emax)           else integer(0)
  idx_ec50_cov <- if (p_ec50 > 0) seq.int(4L + p_lin + p_emax, n_par)                    else integer(0)

  idx <- list(
    e0       = idx_e0,
    emax0    = idx_emax0,
    log_ec50 = idx_log_ec50,
    lin      = idx_lin,
    emax_cov = idx_emax_cov,
    ec50_cov = idx_ec50_cov
  )

  # ---- Parameter names (optimisation scale) ----------------------------------
  theta_names <- c(
    "e0", "emax", "log_ec50",
    if (p_lin  > 0) paste0("linear.",   colnames(Z_lin))  else character(0),
    if (p_emax > 0) paste0("emax.",     colnames(Z_emax)) else character(0),
    if (p_ec50 > 0) paste0("log_ec50.", colnames(Z_ec50)) else character(0)
  )

  # ---- Negative log-likelihood -----------------------------------------------
  nll <- function(theta) {
    emax_i <- theta[idx_emax0] +
      if (p_emax > 0) drop(Z_emax %*% theta[idx_emax_cov]) else 0
    ec50_i <- exp(theta[idx_log_ec50] +
      if (p_ec50 > 0) drop(Z_ec50 %*% theta[idx_ec50_cov]) else 0)
    eta <- theta[idx_e0] + emax_i * x / (ec50_i + x) +
      if (p_lin > 0) drop(Z_lin %*% theta[idx_lin]) else 0
    -sum(
      y       * stats::plogis(eta, log.p = TRUE) +
      (1 - y) * stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
    )
  }

  # ---- Starting values -------------------------------------------------------
  if (is.null(start)) {
    fit_lin <- try(
      suppressWarnings(stats::glm(y ~ x, family = stats::binomial)),
      silent = TRUE
    )
    if (!inherits(fit_lin, "try-error") && all(is.finite(stats::coef(fit_lin)))) {
      b0         <- stats::coef(fit_lin)[[1]]
      b1         <- stats::coef(fit_lin)[[2]]
      emax_start <- b1 * diff(range(x))
      if (abs(emax_start) < 1e-4) emax_start <- 1
      e0_start   <- b0
    } else {
      e0_start   <- stats::qlogis(mean(pmin(pmax(y, 0.01), 0.99)))
      emax_start <- 1
    }
    x_pos      <- unique(x[x > 0])
    ec50_start <- stats::median(x_pos)
    start_full <- c(e0_start, emax_start, log(ec50_start), rep(0, n_par - 3L))
    # Extra EC50 starting points guard against a local optimum.
    ec50_alt   <- stats::quantile(x_pos, c(0.25, 0.75), names = FALSE)

  } else if (length(start) == 3L || length(start) == n_par) {
    start <- as.numeric(start)
    if (!all(is.finite(start)) || start[[3]] <= 0)
      stop("`start` must be finite, with a positive EC50 (third element).")
    start_full <- if (length(start) == 3L)
      c(start[1:2], log(start[[3]]), rep(0, n_par - 3L)) else start
    start_full[idx_log_ec50] <- log(start[[3]])
    ec50_alt <- numeric(0)

  } else {
    stop(sprintf(
      "`start` must have length 3 or %d (= total number of parameters).", n_par
    ))
  }
  names(start_full) <- theta_names

  # ---- Optimise --------------------------------------------------------------
  starts <- c(
    list(start_full),
    lapply(ec50_alt, function(e) replace(start_full, idx_log_ec50, log(e)))
  )
  runs <- lapply(starts, function(s) tryCatch(
    stats::optim(par = s, fn = nll, method = "BFGS", control = control, ...),
    error = function(e) NULL
  ))
  runs <- Filter(function(r) !is.null(r) && is.finite(r$value), runs)
  if (length(runs) == 0L) stop("Optimisation failed from every starting value.")
  opt <- runs[[which.min(vapply(runs, `[[`, numeric(1L), "value"))]]

  # ---- Convergence and identifiability ---------------------------------------
  # The optimiser's status alone isn't enough: also require a positive-definite,
  # well-conditioned Hessian (otherwise some parameter isn't identified) and a
  # negligible Newton decrement, 0.5 * g' H^-1 g, i.e. the log-likelihood the
  # optimiser could still gain. BFGS often stops a little short, so finish
  # with Newton steps first.
  pol       <- .polish_optimum(nll, setNames(opt$par, theta_names))
  th        <- pol$par
  H_inv     <- pol$H_inv
  decrement <- pol$decrement
  nll_hat   <- nll(th)
  ec50_est  <- exp(th[idx_log_ec50])
  converged <- opt$convergence == 0 && !is.null(H_inv) && decrement < 1e-4
  if (opt$convergence != 0) {
    warning("Optimisation did not converge. Code: ", opt$convergence)
  } else if (is.null(H_inv)) {
    .emax_warning(paste0(
      "The Hessian at the optimum is not positive definite or is nearly ",
      "singular, so some parameters are not identified. Standard errors are ",
      "unavailable and `converged` is FALSE."
    ), "bootlace_emax_identification_warning")
  } else if (!converged) {
    .emax_warning(paste0(
      "The optimiser stopped before reaching the optimum; try a smaller ",
      "`control$reltol`."
    ), "bootlace_emax_identification_warning")
  }

  # ---- Is the EC50 estimate supported by the observed exposures? -------------
  ec50_msg <- .ec50_range_message(ec50_est, x)
  if (!is.null(ec50_msg))
    .emax_warning(
      paste0("fit_logistic_emax(): ", ec50_msg),
      "bootlace_emax_ec50_warning"
    )

  # ---- Variance-covariance ---------------------------------------------------
  vcov_theta <- if (hessian) H_inv else NULL
  if (!is.null(vcov_theta)) dimnames(vcov_theta) <- list(theta_names, theta_names)

  # ---- Standard errors (reporting scale) ------------------------------------
  if (!is.null(vcov_theta)) {
    se_th   <- sqrt(diag(vcov_theta))
    std_err <- c(
      se_th[idx_e0], se_th[idx_emax0],
      ec50_est * se_th[idx_log_ec50],           # delta method for ec50
      if (p_lin  > 0) se_th[idx_lin]      else NULL,
      if (p_emax > 0) se_th[idx_emax_cov] else NULL,
      if (p_ec50 > 0) se_th[idx_ec50_cov] else NULL
    )
  } else {
    std_err <- rep(NA_real_, n_par)
  }

  # ---- Coefficients (reporting scale) ----------------------------------------
  coef_names <- c(
    "e0", "emax", "ec50",
    if (p_lin  > 0) paste0("linear.",   colnames(Z_lin))  else character(0),
    if (p_emax > 0) paste0("emax.",     colnames(Z_emax)) else character(0),
    if (p_ec50 > 0) paste0("log_ec50.", colnames(Z_ec50)) else character(0)
  )
  coefs <- setNames(
    c(
      th[idx_e0], th[idx_emax0], ec50_est,
      if (p_lin  > 0) th[idx_lin]      else NULL,
      if (p_emax > 0) th[idx_emax_cov] else NULL,
      if (p_ec50 > 0) th[idx_ec50_cov] else NULL
    ),
    coef_names
  )
  names(std_err) <- coef_names

  # ---- Fitted values ---------------------------------------------------------
  emax_i <- th[idx_emax0] +
    if (p_emax > 0) drop(Z_emax %*% th[idx_emax_cov]) else 0
  ec50_i <- exp(th[idx_log_ec50] +
    if (p_ec50 > 0) drop(Z_ec50 %*% th[idx_ec50_cov]) else 0)
  eta    <- th[idx_e0] + emax_i * x / (ec50_i + x) +
    if (p_lin > 0) drop(Z_lin %*% th[idx_lin]) else 0
  fitted <- as.vector(stats::plogis(eta))

  out <- structure(
    list(
      coefficients      = coefs,
      std_err           = std_err,
      vcov_theta        = vcov_theta,
      theta_hat         = th,
      fitted.values     = fitted,
      linear.predictors = as.vector(eta),
      residuals         = y - fitted,
      loglik            = -nll_hat,
      AIC               = 2 * n_par + 2 * nll_hat,
      converged         = converged,
      n_par             = n_par,
      x                 = x,
      y                 = y,
      Z_lin             = Z_lin,
      Z_emax            = Z_emax,
      Z_ec50            = Z_ec50,
      idx               = idx,
      mm_info_lin       = mm_info_lin,
      mm_info_emax      = mm_info_emax,
      mm_info_ec50      = mm_info_ec50,
      call              = match.call()
    ),
    class = "logistic_emax"
  )
  out$sample_size <- .sample_size(out)
  .warn_sample_size(out$sample_size, "fit_logistic_emax")
  out
}


#' Print a `logistic_emax` object
#'
#' Displays parameter estimates and standard errors, grouped by parameter type.
#'
#' @param x A `logistic_emax` object returned by [fit_logistic_emax()].
#' @param digits Number of decimal places. Default `4`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [fit_logistic_emax()], [predict.logistic_emax()]
#'
#' @export
print.logistic_emax <- function(x, digits = 4, ...) {
  fmt <- function(nms, strip_prefix = NULL) {
    tab <- cbind(Estimate = x$coefficients[nms], `Std. Error` = x$std_err[nms])
    if (!is.null(strip_prefix))
      rownames(tab) <- sub(paste0("^", strip_prefix, "\\."), "", nms)
    round(tab, digits)
  }

  cat("Logistic Emax Model\n\n")
  cat("Core parameters:\n")
  print(fmt(c("e0", "emax", "ec50")), ...)

  lin_nm  <- grep("^linear\\.",   names(x$coefficients), value = TRUE)
  emax_nm <- grep("^emax\\.",     names(x$coefficients), value = TRUE)
  ec50_nm <- grep("^log_ec50\\.", names(x$coefficients), value = TRUE)

  if (length(lin_nm) > 0) {
    cat("\nLinear covariates (additive on link scale):\n")
    print(fmt(lin_nm,  "linear"), ...)
  }
  if (length(emax_nm) > 0) {
    cat("\nEmax covariates (additive effect on Emax):\n")
    print(fmt(emax_nm, "emax"), ...)
  }
  if (length(ec50_nm) > 0) {
    cat("\nEC50 covariates (additive effect on log EC50):\n")
    print(fmt(ec50_nm, "log_ec50"), ...)
  }

  if (!is.null(x$sample_size))
    cat("\nEvents:", .events_line(x$sample_size), "\n")

  cat(
    "\nLog-likelihood:", round(x$loglik, 2),
    "  AIC:", round(x$AIC, 2),
    "  Converged:", x$converged, "\n"
  )
  invisible(x)
}


#' Predict from a `logistic_emax` object
#'
#' Generates predictions on the link, response, or confidence-interval scale.
#' For models with covariates, new covariate matrices can be supplied via
#' `newlinear`, `newemax`, and `newec50`; when these are `NULL` and a new `x`
#' vector is given, covariates are fixed at their training-data column means.
#'
#' Confidence intervals (`type = "ci"`) use the delta method on the link scale,
#' propagating uncertainty through all model parameters including covariate
#' effects.
#'
#' @param object A `logistic_emax` object.
#' @param newdata Optional numeric vector of x values at which to predict.
#'   `NULL` (default) predicts at the original observed values.
#' @param newlinear Optional matrix of linear covariates at the prediction
#'   points. Must have `length(newdata)` rows and the same number of columns as
#'   the training `linear_covs`. When `NULL` and `newdata` is supplied, training
#'   column means are used (with a message).
#' @param newemax Optional matrix of Emax covariates at the prediction points.
#'   Same convention as `newlinear`.
#' @param newec50 Optional matrix of EC50 covariates at the prediction points.
#'   Same convention as `newlinear`.
#' @param type One of `"link"`, `"response"`, or `"ci"`. Default `"link"`.
#' @param level Confidence level for `type = "ci"`. Default `0.95`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return For `type = "link"` or `"response"`, a numeric vector. For
#'   `type = "ci"`, a data frame with columns `x`, `fit`, `lwr`, `upr`.
#'
#' @seealso [fit_logistic_emax()]
#'
#' @examples
#' set.seed(1)
#' dose <- sort(rep(c(0, 5, 10, 25, 50), length.out = 150))
#' eta  <- qlogis(0.10) + 2 * dose / (10 + dose)
#' y    <- rbinom(150, 1, plogis(eta))
#' fit  <- fit_logistic_emax(y, dose)
#'
#' grid <- seq(0, 50, length.out = 100)
#' head(predict(fit, newdata = grid, type = "ci"))
#'
#' @export
predict.logistic_emax <- function(
  object,
  newdata   = NULL,
  newlinear = NULL,
  newemax   = NULL,
  newec50   = NULL,
  type      = c("link", "response", "ci"),
  level     = 0.95,
  ...
) {
  type             <- match.arg(type)
  newdata_provided <- !is.null(newdata)
  x_new            <- if (newdata_provided) as.numeric(newdata) else object$x
  n_new            <- length(x_new)

  # ---- Resolve covariate matrices at prediction points -----------------------
  expand_new_cov <- function(new_z, mm_info) {
    if (is.null(mm_info)) return(as.matrix(new_z))
    df  <- if (is.data.frame(new_z)) new_z else
      data.frame(.var = if (is.character(new_z)) as.factor(new_z) else new_z)
    df[]  <- lapply(df, function(v) if (is.character(v)) as.factor(v) else v)
    mm  <- stats::model.matrix(
      mm_info$terms, data = df,
      xlev          = mm_info$xlevels,
      contrasts.arg = mm_info$contrasts
    )
    mm  <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm), dimnames = dimnames(mm))
  }

  resolve_cov <- function(new_z, train_z, mm_info, label) {
    p <- ncol(train_z)
    if (p == 0L) return(matrix(0, nrow = n_new, ncol = 0L))
    if (!is.null(new_z)) {
      m <- expand_new_cov(new_z, mm_info)
      if (nrow(m) != n_new)
        stop(label, " must have ", n_new, " rows.")
      if (ncol(m) != p)
        stop(label, " must have ", p, " column(s).")
      return(m)
    }
    if (newdata_provided) {
      message(label, " not supplied; fixing at training column means.")
      return(matrix(colMeans(train_z), nrow = n_new, ncol = p, byrow = TRUE))
    }
    train_z
  }

  Z_lin_new  <- resolve_cov(newlinear, object$Z_lin,  object$mm_info_lin,  "newlinear")
  Z_emax_new <- resolve_cov(newemax,   object$Z_emax, object$mm_info_emax, "newemax")
  Z_ec50_new <- resolve_cov(newec50,   object$Z_ec50, object$mm_info_ec50, "newec50")

  p_lin  <- ncol(Z_lin_new)
  p_emax <- ncol(Z_emax_new)
  p_ec50 <- ncol(Z_ec50_new)

  # ---- Linear predictor ------------------------------------------------------
  idx <- object$idx
  th  <- object$theta_hat

  emax_i <- th[idx$emax0] +
    if (p_emax > 0) drop(Z_emax_new %*% th[idx$emax_cov]) else 0
  ec50_i <- exp(th[idx$log_ec50] +
    if (p_ec50 > 0) drop(Z_ec50_new %*% th[idx$ec50_cov]) else 0)
  eta    <- th[idx$e0] + emax_i * x_new / (ec50_i + x_new) +
    if (p_lin > 0) drop(Z_lin_new %*% th[idx$lin]) else 0

  if (type == "link")     return(as.vector(eta))
  if (type == "response") return(as.vector(stats::plogis(eta)))

  # ---- CI via delta method ---------------------------------------------------
  V <- object$vcov_theta
  if (is.null(V)) stop("vcov_theta not available; refit with hessian = TRUE.")

  n_par <- object$n_par
  G     <- matrix(0, nrow = n_new, ncol = n_par)
  denom <- ec50_i + x_new

  G[, idx$e0]       <- 1
  G[, idx$emax0]    <- x_new / denom
  G[, idx$log_ec50] <- -emax_i * x_new * ec50_i / denom^2

  # sweep(M, 1, v, "*") multiplies row i of M by v[i]
  if (p_lin  > 0) G[, idx$lin]      <- Z_lin_new
  if (p_emax > 0) G[, idx$emax_cov] <- sweep(Z_emax_new, 1, x_new / denom, `*`)
  if (p_ec50 > 0) G[, idx$ec50_cov] <- sweep(
    Z_ec50_new, 1, -emax_i * x_new * ec50_i / denom^2, `*`
  )

  var_eta <- rowSums((G %*% V) * G)
  se_eta  <- sqrt(pmax(0, var_eta))
  z       <- stats::qnorm(1 - (1 - level) / 2)

  data.frame(
    x   = x_new,
    fit = stats::plogis(eta),
    lwr = stats::plogis(eta - z * se_eta),
    upr = stats::plogis(eta + z * se_eta)
  )
}


# ---- Internal helpers ---------------------------------------------------------

# Inverse of a Hessian, or NULL when it is not finite, not positive definite,
# or too ill-conditioned to invert reliably (a sign of unidentified
# parameters).
.pd_inverse <- function(H) {
  if (!all(is.finite(H))) return(NULL)
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  if (kappa(H, exact = TRUE) > 1e10) return(NULL)
  chol2inv(R)
}

# Central-difference gradient of `f` at `p`.
.num_grad <- function(f, p, h = 1e-5) {
  vapply(seq_along(p), function(j) {
    e <- numeric(length(p))
    e[j] <- h * max(1, abs(p[[j]]))
    (f(p + e) - f(p - e)) / (2 * e[j])
  }, numeric(1L))
}

# Finish a minimisation of `f` from `p` with Newton steps (with step-halving)
# until the Newton decrement 0.5 * g' H^-1 g is below `tol`. Returns the
# final parameters, the inverse Hessian (NULL when the Hessian isn't positive
# definite) and the decrement.
.polish_optimum <- function(f, p, tol = 1e-4, max_steps = 10L) {
  for (k in 0:max_steps) {
    H_inv <- .pd_inverse(stats::optimHess(p, f))
    if (is.null(H_inv)) return(list(par = p, H_inv = NULL, decrement = NA_real_))
    g         <- .num_grad(f, p)
    step      <- drop(H_inv %*% g)
    decrement <- 0.5 * sum(g * step)
    if (decrement < tol || k == max_steps) break
    f0 <- f(p)
    t  <- 1
    repeat {
      p_new <- p - t * step
      f_new <- f(p_new)
      if (is.finite(f_new) && f_new < f0) break
      t <- t / 2
      if (t < 1e-8) return(list(par = p, H_inv = H_inv, decrement = decrement))
    }
    p <- p_new
  }
  list(par = p, H_inv = H_inv, decrement = decrement)
}

# Raise a warning with an extra condition class, so callers such as
# bootstrap_emax() can tell these apart from fitting failures.
.emax_warning <- function(msg, class) {
  warning(structure(
    class = c(class, "warning", "condition"),
    list(message = msg, call = NULL)
  ))
}

# Explanation of why an EC50 estimate outside the observed non-zero exposure
# range is unreliable, or NULL when it lies inside the range.
.ec50_range_message <- function(ec50, x) {
  lim <- .ec50_support(x)
  hi  <- lim[["hi"]]
  lo  <- lim[["lo"]]
  if (ec50 > hi) {
    sprintf(paste0(
      "the EC50 estimate (%s) is above the highest observed exposure (%s). ",
      "The data don't show where the response levels off, so Emax and EC50 ",
      "aren't determined separately: a larger EC50 with a larger Emax fits ",
      "almost equally well, and their standard errors and confidence ",
      "intervals are unreliable. Interpret only the fitted curve within the ",
      "observed exposure range. Options: a profile-likelihood interval for ",
      "EC50, fixing Emax or EC50 from prior information, a simpler ",
      "exposure-response model (e.g. log-linear), or data at higher exposures."
    ), format(signif(ec50, 3)), format(hi))
  } else if (ec50 < lo) {
    sprintf(paste0(
      "the EC50 estimate (%s) is below a tenth of the lowest non-zero ",
      "exposure (%s). The response is essentially at its maximum already at ",
      "the lowest exposure, so the curve is close to a step and EC50 is set ",
      "by the shape the model assumes rather than by the data; its standard ",
      "error and confidence interval are unreliable. Options: a ",
      "profile-likelihood interval for EC50, fixing EC50 from prior ",
      "information, or data at lower exposures."
    ), format(signif(ec50, 3)), format(10 * lo))
  }
}

# Range of EC50 values the observed exposures can support: from a tenth of
# the lowest non-zero exposure (below that the curve is effectively a step)
# to the highest exposure (above that the plateau is never observed and
# Emax and EC50 trade off).
.EC50_LOW_FRACTION <- 0.1
.ec50_support <- function(x) {
  x_pos <- x[x > 0]
  c(lo = .EC50_LOW_FRACTION * min(x_pos), hi = max(x_pos))
}
