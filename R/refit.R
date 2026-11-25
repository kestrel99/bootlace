# Refitting a model on resampled or reduced data without changing the
# analysis population.
#
# Bootstraps refit from the fitted model's own design (model matrix, response,
# prior weights, offset and, for Cox, strata), so rows excluded by `subset =`
# or missing values never come back and weights and offsets are never lost.
# .analysis_data() gives the matching rows of the user's data for plotting
# and for stepwise refits.

# Rows of the user's data that `fit` was fitted to (after `subset =` and
# missing-value removal), or NULL when that data can't be recovered.
.analysis_data <- function(fit) {
  data <- tryCatch(
    eval(fit$call$data, envir = environment(stats::formula(fit))),
    error = function(e) NULL
  )
  if (!is.data.frame(data)) return(NULL)
  rows <- rownames(stats::model.frame(fit))
  if (!all(rows %in% rownames(data))) return(NULL)
  data[rows, , drop = FALSE]
}

# ---- Binomial GLM -------------------------------------------------------------

.glm_design <- function(fit) {
  X <- stats::model.matrix(fit)
  list(
    n       = nrow(X),
    X       = X,
    y       = fit$y,
    weights = fit$prior.weights,
    offset  = fit$offset,
    family  = fit$family,
    control = fit$control,
    cluster = NULL
  )
}

# Coefficients of the GLM refitted to rows `idx` of its design. Columns that
# are all zero in a resample (e.g. an absent factor level) come back as NA.
.refit_glm <- function(des, idx) {
  f <- stats::glm.fit(
    x       = des$X[idx, , drop = FALSE],
    y       = des$y[idx],
    weights = des$weights[idx],
    offset  = des$offset[idx],
    family  = des$family,
    control = des$control
  )
  stats::setNames(f$coefficients, colnames(des$X))
}

# ---- Cox ------------------------------------------------------------------------

.cox_design <- function(fit) {
  mf <- stats::model.frame(fit)
  tt <- stats::terms(fit)
  y  <- fit$y
  if (is.null(y)) y <- stats::model.response(mf)

  st_vars <- survival::untangle.specials(tt, "strata")$vars
  strata  <- if (length(st_vars) > 0L)
    survival::strata(mf[st_vars], shortlabel = TRUE) else NULL

  cl_vars <- survival::untangle.specials(tt, "cluster")$vars
  cluster <- if (!is.null(mf[["(id)"]])) {
    mf[["(id)"]]
  } else if (!is.null(mf[["(cluster)"]])) {
    mf[["(cluster)"]]
  } else if (length(cl_vars) > 0L) {
    mf[[cl_vars[1L]]]
  } else {
    NULL
  }

  list(
    n       = nrow(y),
    X       = stats::model.matrix(fit),
    y       = y,
    weights = stats::model.weights(mf),
    offset  = stats::model.offset(mf),
    strata  = strata,
    cluster = cluster,
    method  = fit$method,
    type    = attr(y, "type")
  )
}

# Log hazard ratios of the Cox model refitted to rows `idx` of its design.
# With `basehaz = TRUE` (unstratified designs only), returns
# list(coef =, extra =) where `extra` is the refit's uncentred baseline
# cumulative hazard as a two-column matrix (time, hazard).
.refit_cox <- function(des, idx, basehaz = FALSE) {
  env <- list2env(
    list(
      Y = des$y[idx], X = des$X[idx, , drop = FALSE],
      off = des$offset[idx], st = des$strata[idx], W = des$weights[idx],
      strata = survival::strata, offset = stats::offset
    ),
    parent = baseenv()
  )
  rhs <- c(
    "X",
    if (!is.null(des$offset)) "offset(off)",
    if (!is.null(des$strata)) "strata(st)"
  )
  f <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")), env = env)
  # coxph() evaluates `weights` in the formula's environment, where W lives.
  m <- survival::coxph(f, weights = W, ties = des$method)
  b <- stats::setNames(as.numeric(stats::coef(m)), colnames(des$X))
  if (!basehaz || !is.null(des$strata)) return(b)
  h0 <- survival::basehaz(m, centered = FALSE)
  list(coef = b, extra = cbind(time = h0$time, hazard = h0$hazard))
}

# ---- Resampling -------------------------------------------------------------------

# Bootstrap row indices into a design: whole clusters (subjects) when the
# design has a cluster identifier, otherwise single rows.
.resample_rows <- function(des) {
  if (is.null(des$cluster)) return(sample.int(des$n, des$n, replace = TRUE))
  by_cluster <- split(seq_len(des$n), des$cluster)
  drawn <- sample.int(length(by_cluster), length(by_cluster), replace = TRUE)
  unlist(by_cluster[drawn], use.names = FALSE)
}
