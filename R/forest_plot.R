# Internal formatting helpers --------------------------------------------------

# Format to `digits` significant figures; trimws() removes the leading space
# that formatC inserts when right-aligning numeric output.
.fp_fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "-",
         trimws(formatC(signif(x, digits), format = "g", digits = digits)))
}

# Format p-values: "" for NA, "<0.001" for very small values.
.fp_fmt_p <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

# fp_data ----------------------------------------------------------------------

#' Build a tidy data frame for a forest plot
#'
#' Prepares a structured data frame from a fitted binomial `glm` or `coxph`
#' object suitable for passing to [fp_plot()].  The model type is
#' auto-detected from the object class.
#'
#' For categorical predictors one row is produced per factor level, with the
#' first level as the reference category.  The p-value on the header row is
#' derived from the Wald statistic of the first non-reference contrast
#' (appropriate for binary predictors; compute a likelihood-ratio p-value
#' separately for predictors with more than two levels).
#'
#' For continuous predictors each row represents an estimated OR/HR relative
#' to the predictor's median in `data`.
#'
#' Interaction terms (those containing `":"`) are detected automatically.  A
#' continuous x continuous interaction is shown as a single row.  An
#' interaction involving a factor is shown as a header row with one child row
#' per factor level present in the model.  The exponentiated coefficient of an
#' interaction is a ratio of odds ratios (or of hazard ratios), not an OR/HR
#' itself, so these rows are labelled "(ratio of ORs)" or "(ratio of HRs)".
#'
#' Every displayed estimate is a linear combination of the fitted
#' coefficients built from the coding R actually used for each factor in
#' each term (contrasts, or indicator coding when a lower-order term is
#' absent, as for `a` in `y ~ a + a:b`), so the table is the same whatever
#' contrast coding was used (e.g. `contr.treatment` or `contr.sum`). Factor rows compare each level with the first (reference)
#' level. When a term is part of an interaction, its rows are its effect at
#' the reference level of each factor it interacts with (and at 0 for an
#' interacting numeric variable), which is what treatment coding reports.
#' Factor-interaction rows are ratios of ORs/HRs for each combination of
#' non-reference levels.
#'
#' Header-row p-values test the whole term: a joint Wald test that all of
#' those comparisons are zero. For a factor with more than two levels, or a
#' factor interaction with several contrasts, the P-value column shows the
#' test's degrees of freedom, e.g. `0.023 (2 df)`. Single-coefficient terms
#' get the ordinary Wald test. A term with no estimable coefficient (fully
#' aliased) shows "not estimable" and no level rows.
#'
#' @param model   A fitted `coxph` object, or a `glm` with a binomial (or
#'   quasibinomial) family and logit link. Other models are rejected, because
#'   their exponentiated coefficients are not hazard or odds ratios.
#' @param data    Data frame used to fit `model`.
#' @param labels  Named character vector mapping raw variable names to
#'   publication-ready display labels.  Unmatched names fall back to the
#'   raw variable name.  Interaction term names (e.g. `"age:bmi"`) can also
#'   be supplied as keys.
#' @param log_vars Character vector of variable names stored on the log scale.
#'   Their percentile display values are back-transformed with `exp()` so
#'   readers see the original units.
#' @param percs   Numeric vector of percentiles shown for continuous
#'   predictors (default `c(0.05, 0.25, 0.75, 0.95)`).
#' @param ci_level Confidence level for the CI text column (default `0.95`).
#' @param digits  Significant figures used when formatting percentile display
#'   values and the estimate/CI text column (default `3`).
#' @param ci_width Integer; number of spaces in the CI-panel spacer column.
#'   Increase to widen the CI drawing area in [fp_plot()] (default `34`).
#'
#' @return A data frame with one row per displayed entry containing:
#'   `Predictor`, `est`, `lci`, `uci`, `pval`, `p_df` (degrees of freedom of
#'   the test behind `pval`), `is_header`, `is_reference`,
#'   `group`, the model-specific estimate column (`OR (95% CI)` or
#'   `HR (95% CI)`), `P-value`, and `" "` (spacer column).
#'
#' @seealso [fp_plot()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   lung$sex <- factor(lung$sex, levels = 1:2, labels = c("Male", "Female"))
#'   lung_cc <- lung[
#'     complete.cases(lung[, c("time", "status", "age", "sex")]), ]
#'
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex,
#'     data = lung_cc
#'   )
#'   pd <- fp_data(fit, lung_cc,
#'                 labels = c(age = "Age (years)", sex = "Sex"))
#'   head(pd[, c("Predictor", "HR (95% CI)", "P-value")])
#' }
#'
#' @export
fp_data <- function(
  model,
  data,
  labels   = NULL,
  log_vars = NULL,
  percs    = c(0.05, 0.25, 0.75, 0.95),
  ci_level = 0.95,
  digits   = 3L,
  ci_width = 34L
) {
  .check_prob(ci_level, "ci_level")
  .check_range(percs, "percs", lower = 0, upper = 1)

  # exp(coef) is a hazard ratio only for Cox models and an odds ratio only for
  # a logit link; label nothing else as either.
  is_cox   <- inherits(model, "coxph")
  is_logit <- inherits(model, "glm") &&
    model$family$family %in% c("binomial", "quasibinomial") &&
    identical(model$family$link, "logit")
  if (!is_cox && !is_logit)
    stop(
      "`model` must be a coxph fit or a binomial glm with a logit link. ",
      "Exponentiated coefficients of other models are not hazard or odds ",
      "ratios.",
      call. = FALSE
    )
  est_label <- if (is_cox) "HR (95% CI)" else "OR (95% CI)"

  coefs      <- stats::coef(model)
  vcov_mat   <- stats::vcov(model)
  coef_names <- names(coefs)
  z_crit     <- stats::qnorm(1 - (1 - ci_level) / 2)
  terms_     <- attr(stats::terms(model), "term.labels")

  # -- Input validation ------------------------------------------------------
  if (!is.null(labels)) {
    bad_lbl <- setdiff(names(labels), terms_)
    if (length(bad_lbl))
      warning("fp_data(): label key(s) not found in model terms: ",
              paste(bad_lbl, collapse = ", "),
              call. = FALSE)
  }

  if (!is.null(log_vars)) {
    bad_lv <- setdiff(log_vars, terms_)
    if (length(bad_lv))
      warning("fp_data(): log_vars entry/entries not found in model terms: ",
              paste(bad_lv, collapse = ", "),
              call. = FALSE)
    cat_lv <- intersect(log_vars, terms_)
    cat_lv <- cat_lv[vapply(cat_lv, function(v) {
      col <- data[[v]]
      !is.null(col) && (is.factor(col) || is.character(col))
    }, logical(1L))]
    if (length(cat_lv))
      warning("fp_data(): log_vars contains categorical predictor(s) ",
              "(back-transform has no effect): ",
              paste(cat_lv, collapse = ", "),
              call. = FALSE)
  }

  .lbl  <- function(v) {
    if (!is.null(labels) && v %in% names(labels)) labels[[v]] else v
  }
  .ilog <- function(v) !is.null(log_vars) && v %in% log_vars

  # Every displayed effect is an explicit linear combination L'b of the
  # fitted coefficients, built from the fitted contrasts (see
  # .fp_effect_L()), so the table doesn't depend on the contrast coding.
  # Effects of a term are taken at the reference level of any factor it
  # interacts with (and at 0 for interacting numeric variables), which is
  # what treatment coding reports.
  term_cols <- .fp_term_cols(model)
  info      <- .fp_design_info(model, terms_, term_cols)

  # Estimate and SE of L'b; NA when L involves a non-estimable coefficient.
  .lin <- function(L) {
    if (is.null(L) || any(L != 0 & is.na(coefs))) return(c(NA_real_, NA_real_))
    nz <- L != 0
    c(sum(L[nz] * coefs[nz]),
      sqrt(drop(crossprod(L[nz], vcov_mat[nz, nz, drop = FALSE] %*% L[nz]))))
  }
  # Joint Wald test of the rows of L (each a combination); NA (df 0) when not
  # estimable.
  .wald_L <- function(Ls) {
    Ls <- Filter(function(L) !is.null(L) && !any(L != 0 & is.na(coefs)), Ls)
    if (length(Ls) == 0L) return(list(p = NA_real_, df = 0L))
    M  <- do.call(rbind, Ls)
    nz <- colSums(M != 0) > 0
    M  <- M[, nz, drop = FALSE]
    est <- drop(M %*% coefs[nz])
    V   <- M %*% vcov_mat[nz, nz, drop = FALSE] %*% t(M)
    w   <- tryCatch(drop(crossprod(est, solve(V, est))), error = function(e) NA_real_)
    if (!is.finite(w)) return(list(p = NA_real_, df = 0L))
    list(p = stats::pchisq(w, df = nrow(M), lower.tail = FALSE), df = nrow(M))
  }

  # exp() of an interaction coefficient is a ratio of ORs / HRs.
  ratio_lbl <- if (is_cox) " (ratio of HRs)" else " (ratio of ORs)"

  .row <- function(Predictor, est = NA_real_, lci = NA_real_, uci = NA_real_,
                   pval = NA_real_, p_df = NA_integer_, is_header = FALSE,
                   is_reference = FALSE) {
    data.frame(
      Predictor = Predictor, est = est, lci = lci, uci = uci, pval = pval,
      p_df = as.integer(p_df), is_header = is_header,
      is_reference = is_reference, group = grp, stringsAsFactors = FALSE
    )
  }
  .exp_row <- function(Predictor, es) {
    .row(Predictor, est = exp(es[1]), lci = exp(es[1] - z_crit * es[2]),
         uci = exp(es[1] + z_crit * es[2]))
  }

  rows <- list()
  grp  <- 0L

  for (term in terms_) {
    grp <- grp + 1L
    if (length(term_cols[[term]]) == 0L) next
    vars_t <- info$vars[[term]]
    facs   <- vars_t[vars_t %in% info$factors]
    is_interaction <- length(vars_t) > 1L

    if (length(facs) > 0L) {
      # -- Factor main effect or interaction involving a factor ---------------
      cmps <- .fp_comparisons(info, facs)
      Ls   <- lapply(cmps$idx, function(ix) .fp_effect_L(info, term, ix, coef_names))
      wt   <- .wald_L(Ls)
      hdr  <- if (is_interaction) paste0(.lbl(term), ratio_lbl) else .lbl(term)
      rows <- c(rows, list(.row(hdr, pval = wt$p, p_df = wt$df,
                                is_header = TRUE)))
      if (!is_interaction)
        rows <- c(rows, list(.row(
          paste0("  ", info$levels[[term]][[1L]], " (reference)"),
          is_reference = TRUE
        )))
      for (k in seq_along(Ls)) {
        es <- .lin(Ls[[k]])
        if (is.na(es[1])) next
        rows <- c(rows, list(.exp_row(paste0("  ", cmps$label[[k]]), es)))
      }

    } else if (is_interaction) {
      # -- Numeric x numeric interaction: one row -----------------------------
      es <- .lin(.fp_effect_L(info, term, list(), coef_names))
      rows <- c(rows, list(.row(
        paste0(.lbl(term), ratio_lbl),
        est = exp(es[1]), lci = exp(es[1] - z_crit * es[2]),
        uci = exp(es[1] + z_crit * es[2]),
        pval = 2 * stats::pnorm(-abs(es[1] / es[2])), p_df = 1L
      )))

    } else {
      # -- Continuous main effect: header + percentile rows -------------------
      col <- data[[term]]
      if (!is.numeric(col)) next     # e.g. a logical predictor: not supported
      es  <- .lin(.fp_effect_L(info, term, list(), coef_names))
      b   <- es[1]
      s   <- es[2]
      p   <- if (is.na(b)) NA_real_ else 2 * stats::pnorm(-abs(b / s))
      rows <- c(rows, list(.row(.lbl(term), pval = p,
                                p_df = if (is.na(b)) 0L else 1L,
                                is_header = TRUE)))
      if (is.na(b)) next
      med <- stats::median(col, na.rm = TRUE)
      for (perc in percs) {
        pv      <- stats::quantile(col, probs = perc, na.rm = TRUE)
        d       <- pv - med
        dv      <- if (.ilog(term)) exp(pv) else pv
        pct_lbl <- sprintf("%dth percentile", as.integer(perc * 100))
        r <- .exp_row(paste0("  ", pct_lbl, " (", .fp_fmt(dv, digits), ")"),
                      c(b * d, s * abs(d)))
        r$is_reference <- abs(d) < 1e-10
        rows <- c(rows, list(r))
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  ci_str <- ifelse(
    is.na(out$est),
    "",
    paste0(
      .fp_fmt(out$est, digits), " (",
      .fp_fmt(out$lci, digits), ", ",
      .fp_fmt(out$uci, digits), ")"
    )
  )

  out[[est_label]] <- ci_str
  out[["P-value"]] <- paste0(
    .fp_fmt_p(out$pval),
    ifelse(!is.na(out$p_df) & out$p_df > 1L, paste0(" (", out$p_df, " df)"), "")
  )
  # Header tests with no estimable coefficient
  not_est <- out$is_header & is.na(out$pval) & !is.na(out$p_df) & out$p_df == 0L
  out[["P-value"]][not_est] <- "not estimable"
  out[[" "]]       <- strrep(" ", ci_width)

  out
}


# fp_plot ----------------------------------------------------------------------

#' Build and render a forest plot from a fitted model
#'
#' Calls [fp_data()] to prepare the plot data frame and then renders a
#' `forestploter` forest plot in one step.  All [fp_data()] arguments are
#' forwarded directly; call [fp_data()] separately first only if you need to
#' inspect or modify the intermediate data frame before plotting.
#'
#' The model type (logistic or Cox) is auto-detected from the object class.
#'
#' The four display columns are:
#' 1. `Predictor`  -- row labels
#' 2. `OR (95% CI)` or `HR (95% CI)`  -- formatted estimate string
#' 3. `" "`  -- spacer column where the CI glyph is drawn
#' 4. `P-value`
#'
#' Header rows (one per predictor) receive a shaded background and bold label.
#' P-values below 0.05 are coloured red.
#'
#' @param model   A fitted `coxph` object, or a logit-link binomial `glm`; see
#'   [fp_data()].
#' @param data    Data frame used to fit `model`.
#' @param labels  Named character vector; display labels for predictor names
#'   (passed to [fp_data()]).
#' @param log_vars Character vector of variable names stored on the log scale
#'   (passed to [fp_data()]).
#' @param percs   Percentiles shown for continuous predictors
#'   (passed to [fp_data()], default `c(0.05, 0.25, 0.75, 0.95)`).
#' @param ci_level Confidence level (passed to [fp_data()], default `0.95`).
#' @param digits  Significant figures for formatted text
#'   (passed to [fp_data()], default `3`).
#' @param ci_width Spacer column width (passed to [fp_data()], default `34`).
#' @param xlim       Length-2 numeric; x-axis limits on the original
#'   (non-transformed) scale (default `c(0.1, 10)`).
#' @param vert_line  Numeric vector of additional vertical reference lines
#'   (default `c(0.8, 1.25)`).  Set to `NULL` to suppress.
#' @param ci_col     Colour string for CI points and bars (default dark red).
#' @param base_size  Base font size (default `10`).
#' @param x_trans    Axis transformation: `"log"` (default), `"log2"`,
#'   `"log10"`, or `"none"`.
#' @param ticks_at   Numeric vector of custom tick positions on the original
#'   x-axis scale.  `NULL` (default) uses automatic ticks.
#'
#' @return A `gtable` object of class `"forestplot"`.  Print it directly or
#'   draw it with `grid::grid.draw()` inside a graphics device.
#'
#' @seealso [fp_data()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   lung$sex <- factor(lung$sex, levels = 1:2, labels = c("Male", "Female"))
#'   lung_cc <- lung[
#'     complete.cases(lung[, c("time", "status", "age", "sex")]), ]
#'
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex,
#'     data = lung_cc
#'   )
#'   fp_plot(fit, lung_cc,
#'           labels = c(age = "Age (years)", sex = "Sex"))
#' }
#'
#' @export
fp_plot <- function(
  model,
  data,
  labels    = NULL,
  log_vars  = NULL,
  percs     = c(0.05, 0.25, 0.75, 0.95),
  ci_level  = 0.95,
  digits    = 3L,
  ci_width  = 34L,
  xlim      = c(0.1, 10),
  vert_line = c(0.8, 1.25),
  ci_col    = "#C00000",
  base_size = 10,
  x_trans   = "log",
  ticks_at  = NULL
) {
  plot_data <- fp_data(
    model    = model,
    data     = data,
    labels   = labels,
    log_vars = log_vars,
    percs    = percs,
    ci_level = ci_level,
    digits   = digits,
    ci_width = ci_width
  )

  is_cox   <- "HR (95% CI)" %in% names(plot_data)
  est_col  <- if (is_cox) "HR (95% CI)" else "OR (95% CI)"
  arr_labs <- if (is_cox)
    c("Lower hazard", "Greater hazard")
  else
    c("Lower likelihood", "Greater likelihood")

  disp_data <- plot_data[, c("Predictor", est_col, " ", "P-value")]

  tm <- forestploter::forest_theme(
    base_size        = base_size,
    summary_col      = "black",
    arrow_label_just = "end",
    arrow_type       = "closed",
    ci_pch           = 16,
    ci_col           = ci_col,
    ci_fill          = ci_col,
    ci_alpha         = 1,
    ci_lty           = 1,
    ci_lwd           = 1.5,
    ci_Theight       = 0.2
  )

  fp <- forestploter::forest(
    disp_data,
    est       = plot_data$est,
    lower     = plot_data$lci,
    upper     = plot_data$uci,
    ci_column = 3L,
    ref_line  = 1,
    vert_line = vert_line,
    sizes     = 1,
    arrow_lab = arr_labs,
    x_trans   = x_trans,
    xlim      = xlim,
    ticks_at  = ticks_at,
    theme     = tm
  )

  hrows <- which(plot_data$is_header)
  if (length(hrows) > 0L) {
    fp <- forestploter::edit_plot(
      fp, col = seq_len(ncol(disp_data)), row = hrows,
      which = "background", gp = grid::gpar(fill = "ivory2")
    )
    fp <- forestploter::edit_plot(
      fp, col = 1L, row = hrows,
      gp = grid::gpar(fontface = "bold")
    )
  }

  sig_rows <- which(
    plot_data$is_header & !is.na(plot_data$pval) & plot_data$pval < 0.05
  )
  if (length(sig_rows) > 0L) {
    fp <- forestploter::edit_plot(fp, col = 4L, row = sig_rows,
                                  gp = grid::gpar(col = "red"))
  }

  forestploter::add_border(fp, part = "header", row = 1L,
                           gp = grid::gpar(lwd = 1))
}


# Coefficient names belonging to each model term, from the fit's own
# term-to-column map (coxph: `assign`; glm: the model matrix's `assign`
# attribute), so nothing depends on how contrast columns are named.
.fp_term_cols <- function(model) {
  tl <- attr(stats::terms(model), "term.labels")
  if (inherits(model, "coxph")) {
    cn  <- names(stats::coef(model))
    out <- lapply(model$assign, function(i) cn[i])
    return(out[intersect(tl, names(out))])
  }
  mm   <- stats::model.matrix(model)
  asg  <- attr(mm, "assign")
  keep <- asg > 0
  split(colnames(mm)[keep], factor(tl[asg[keep]], levels = tl))
}

# Contrast matrix of a factor term as fitted, one row per level.
.fp_contrast_matrix <- function(model, term, levs) {
  ctr <- model$contrasts[[term]]
  if (is.null(ctr)) ctr <- "contr.treatment"
  C <- if (is.character(ctr)) do.call(ctr, list(levs)) else as.matrix(ctr)
  rownames(C) <- levs
  C
}


# Design information for building effects: the variables of each term, which
# are factors, their levels and fitted contrast matrices, and how each term's
# columns are coded (fitted contrasts, or full indicators when a lower-order
# term is missing).
.fp_design_info <- function(model, terms_, term_cols) {
  xlev  <- model$xlevels
  vars  <- stats::setNames(lapply(terms_, function(t) strsplit(t, ":", fixed = TRUE)[[1L]]),
                           terms_)
  facs  <- names(xlev)
  C     <- stats::setNames(lapply(facs, function(v) .fp_contrast_matrix(model, v, xlev[[v]])),
                           facs)
  levels <- stats::setNames(lapply(terms_, function(t) xlev[[t]]), terms_)
  # Per term: for each factor in it, the coding matrix (rows = levels) that R
  # actually used for that factor in that term. The terms object's "factors"
  # attribute records it: 1 = contrasts, 2 = indicator (dummy) coding, which
  # R uses when a lower-order term is absent (e.g. `a` in `y ~ a + a:b`).
  fa <- attr(stats::terms(model), "factors")
  coding <- stats::setNames(lapply(terms_, function(t) {
    fv <- vars[[t]][vars[[t]] %in% facs]
    if (length(fv) == 0L) return(list())
    cod <- stats::setNames(lapply(fv, function(v) {
      if (!is.null(fa) && v %in% rownames(fa) && t %in% colnames(fa) &&
          fa[v, t] == 2L) {
        k <- length(xlev[[v]])
        m <- diag(k)
        dimnames(m) <- list(xlev[[v]], xlev[[v]])
        m
      } else {
        C[[v]]
      }
    }), fv)
    # Numeric variables contribute one column each, so the factor codings
    # must account for all of the term's fitted columns.
    if (prod(vapply(cod, ncol, numeric(1L))) != length(term_cols[[t]]))
      return(NULL)
    cod
  }), terms_)
  list(vars = vars, factors = facs, levels = levels, xlevels = xlev,
       coding = coding, term_cols = term_cols, terms = terms_)
}

# Non-reference level combinations of the factors `facs` (first varies
# fastest, matching model-matrix column order), with display labels.
.fp_comparisons <- function(info, facs) {
  grid <- expand.grid(lapply(facs, function(v) seq.int(2L, length(info$xlevels[[v]]))),
                      KEEP.OUT.ATTRS = FALSE)
  names(grid) <- facs
  idx <- lapply(seq_len(nrow(grid)), function(r) as.list(grid[r, , drop = FALSE]))
  label <- vapply(idx, function(ix) {
    paste(vapply(facs, function(v) info$xlevels[[v]][[ix[[v]]]], character(1L)),
          collapse = ":")
  }, character(1L))
  list(idx = idx, label = label)
}

# Coefficient vector L (named by `coef_names`) for the effect of `term`:
# level ix[[v]] vs the reference for each factor v of the term, per unit of
# each numeric variable, evaluated at the reference level of every other
# factor (and 0 for every other numeric variable) in higher-order terms that
# contain it. NULL when a needed term's coding isn't recognised.
.fp_effect_L <- function(info, term, ix, coef_names) {
  L  <- stats::setNames(numeric(length(coef_names)), coef_names)
  tv <- info$vars[[term]]
  for (h in info$terms) {
    hv <- info$vars[[h]]
    if (!all(tv %in% hv)) next
    extra <- setdiff(hv, tv)
    if (any(!extra %in% info$factors)) next     # a numeric extra is set to 0
    coding <- info$coding[[h]]
    if (is.null(coding)) return(NULL)
    vecs <- lapply(hv, function(v) {
      if (!v %in% info$factors) return(1)
      Cv <- coding[[v]]
      if (v %in% tv) Cv[ix[[v]], ] - Cv[1L, ] else Cv[1L, ]
    })
    Lh <- Reduce(function(acc, v) kronecker(v, acc), vecs)
    cols <- info$term_cols[[h]]
    if (length(Lh) != length(cols)) return(NULL)
    L[cols] <- L[cols] + Lh
  }
  L
}
