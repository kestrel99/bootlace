#' Visual predictive check for a Cox proportional hazards model
#'
#' Simulates an event time for every subject in the fitting data from the
#' fitted model: the Breslow baseline cumulative hazard scaled by the subject's
#' hazard ratio, inverted at a uniform draw. The simulated data are censored
#' as the observed data were, and the observed Kaplan-Meier curve is compared
#' with the median and prediction interval of the simulated Kaplan-Meier
#' curves. With `bootstrap = TRUE`, each replicate refits the model to a
#' bootstrap sample of subjects, so the interval also reflects uncertainty in
#' the coefficients and the baseline hazard.
#'
#' @section Censoring:
#' A censored subject keeps their observed censoring time. A subject who had
#' the event never had theirs observed, so it is either `t_max`
#' (`censoring = "admin"`, right when all censoring is administrative) or a
#' draw from the reverse Kaplan-Meier censoring distribution, conditional on
#' exceeding their event time (`"reverseKM"`, for dropout). A simulated event
#' beyond the last observed event time can't be placed by a Cox model, so it
#' is censored too.
#'
#' @param fit A [survival::coxph()] fit to right-censored data, without
#'   `strata()` terms.
#' @param data The data `fit` was fitted to, in the same row order. Rows the
#'   fit dropped for missing values are dropped here too. The remaining rows
#'   must reproduce the fit's response and linear predictor row by row;
#'   reordered or altered data are an error, since the simulated and
#'   observed curves would otherwise pair the wrong records.
#' @param by Name of a discrete column in `data` to stratify the VPC by (e.g.
#'   an exposure quartile). `NULL` for a single panel.
#' @param n_sim Number of simulated replicates.
#' @param pi_probs Lower and upper quantiles of the prediction interval.
#' @param t_max End of follow-up: the latest possible censoring time. Defaults
#'   to the last observed time. If earlier, the observed data are censored at
#'   `t_max` as well.
#' @param n_grid Number of times the simulated curves are evaluated at.
#' @param seed Random seed, or `NULL` to leave the RNG state alone.
#' @param bootstrap Refit the model to a bootstrap sample in each replicate.
#'   When the fit has a subject or cluster identifier (`id =`, `cluster =` or
#'   a `cluster()` term), whole clusters are resampled. Replicates whose refit
#'   or prediction fails (e.g. a sample missing a rare factor level) are
#'   dropped with a warning.
#' @param censoring `"admin"` or `"reverseKM"`; see the Censoring section.
#' @param overlay Draw the `by` groups on one set of axes, one colour per
#'   group, instead of one panel per group. Up to 8 groups.
#'
#' @return A list with:
#'   \describe{
#'     \item{`sim`}{Simulated median and prediction interval (`lo`, `median`,
#'       `hi`) by `group` and `time`.}
#'     \item{`obs`}{Observed Kaplan-Meier curve by `group`, with a `censored`
#'       flag for times with censoring but no events.}
#'     \item{`zph`}{The [survival::cox.zph()] test of proportional hazards, or
#'       `NULL` if it couldn't be run.}
#'     \item{`sample_size`}{Sample-size diagnostics for `fit`; see
#'       **Sample size checks**.}
#'     \item{`plot`}{A `ggplot` object. Its subtitle gives the Schoenfeld
#'       residual p-values: global, and for `by` when it is a model term.}
#'   }
#'
#' @seealso [bootstrap_cox()], [plot_cox_survival()]
#'
#' @examples
#' lung <- survival::lung
#' lung$status <- lung$status - 1L
#' fit <- survival::coxph(
#'   survival::Surv(time, status) ~ age + sex + ph.ecog,
#'   data = lung
#' )
#' vpc <- vpc_cox(fit, lung, by = "sex", n_sim = 100, seed = 1)
#' vpc$plot
#'
#' @template sample-size
#' @export
vpc_cox <- function(
  fit,
  data,
  by = NULL,
  n_sim = 500,
  pi_probs = c(0.05, 0.95),
  t_max = NULL,
  n_grid = 200,
  seed = NULL,
  bootstrap = FALSE,
  censoring = c("admin", "reverseKM"),
  overlay = FALSE
) {
  censoring <- match.arg(censoring)
  n_sim  <- .check_count(n_sim, "n_sim")
  n_grid <- .check_count(n_grid, "n_grid", min = 2L)
  .check_range(pi_probs, "pi_probs", lower = 0, upper = 1)
  if (length(pi_probs) != 2L || !(0 < pi_probs[1] && pi_probs[1] < pi_probs[2] &&
                                  pi_probs[2] < 1)) {
    stop("`pi_probs` must be two increasing numbers strictly between 0 and 1",
         call. = FALSE)
  }
  if (!is.null(t_max) && (length(t_max) != 1L || is.na(t_max) ||
                          !is.finite(t_max) || t_max <= 0)) {
    stop("`t_max` must be a single positive number or NULL", call. = FALSE)
  }
  .check_seed(seed)
  if (!inherits(fit, "coxph")) {
    stop("`fit` must be a coxph() fit", call. = FALSE)
  }
  if (is.null(fit$y)) {
    stop("refit with coxph(..., y = TRUE)", call. = FALSE)
  }
  if (attr(fit$y, "type") != "right") {
    stop(
      "only right-censored data are supported, not (start, stop] data",
      call. = FALSE
    )
  }
  if (!is.null(attr(fit$terms, "specials")$strata)) {
    stop("models with strata() aren't supported", call. = FALSE)
  }
  if (!is.null(fit$na.action)) {
    data <- data[-fit$na.action, , drop = FALSE]
  }
  if (nrow(data) != nrow(fit$y)) {
    stop("`data` doesn't match the data `fit` was fitted to", call. = FALSE)
  }
  # Predictions follow the order of `data` while the observed times follow the
  # fit's, so the rows must be the same records in the same order.
  if (!.vpc_data_matches(fit, data)) {
    stop(
      "`data` must hold the rows `fit` was fitted to, in the same order and ",
      "with the same values",
      call. = FALSE
    )
  }
  if (!is.null(by) && !by %in% names(data)) {
    stop("`data` has no column ", by, call. = FALSE)
  }
  if (
    overlay &&
      !is.null(by) &&
      length(unique(data[[by]])) > length(.vpc_group_colours)
  ) {
    stop(
      "`overlay` supports up to ",
      length(.vpc_group_colours),
      " groups; use panels instead",
      call. = FALSE
    )
  }
  sample_size <- .sample_size(fit)
  .warn_sample_size(sample_size, "vpc_cox")
  if (!is.null(seed)) {
    set.seed(seed)
  }

  time <- fit$y[, "time"]
  status <- fit$y[, "status"]
  n <- length(time)
  if (is.null(t_max)) {
    t_max <- max(time)
  }
  cens <- pmin(time, t_max)
  draw_censoring <- if (censoring == "admin") {
    function() rep(t_max, sum(status == 1))
  } else {
    .reverse_km_sampler(time, status, t_max)
  }

  group <- if (is.null(by)) rep("All", n) else data[[by]]
  group <- factor(group)
  levels(group) <- sprintf("%s (n = %d)", levels(group), tabulate(group))

  grid <- seq(0, t_max, length.out = n_grid)

  # Bootstrap whole subjects/clusters when the fit has an id or cluster.
  boot_design <- list(n = n, cluster = .cox_design(fit)$cluster)

  sim_once <- function() {
    # Uncentred baseline cumulative hazard and linear predictor, so that
    # S_i(t) = exp(-H0(t) * risk_i), for the original subjects.
    f <- fit
    if (bootstrap) {
      boot <- data[.resample_rows(boot_design), , drop = FALSE]
      # model = TRUE keeps the model frame; without it predict() can
      # re-evaluate `data = boot` in the formula's environment and not find it
      # (x = TRUE alone isn't enough, e.g. when the fit has weights).
      f <- tryCatch(
        suppressWarnings(stats::update(fit, . ~ ., data = boot, model = TRUE)),
        error = function(e) paste("refit:", conditionMessage(e))
      )
      if (is.character(f)) return(f)
    }
    risk <- tryCatch(
      exp(stats::predict(f, newdata = data, type = "lp", reference = "zero")),
      error = function(e) paste("predict:", conditionMessage(e))
    )
    if (is.character(risk)) {
      return(risk)
    }
    bh <- survival::basehaz(f, centered = FALSE)

    target <- -log(stats::runif(n)) / risk
    # First time at which H0 reaches the target; beyond the last one: no event.
    idx <- findInterval(target, bh$hazard, left.open = TRUE) + 1
    t_event <- rep(Inf, n)
    t_event[idx <= nrow(bh)] <- bh$time[idx[idx <= nrow(bh)]]
    cens[status == 1] <- draw_censoring()
    .km_at_grid(
      pmin(t_event, cens),
      as.integer(t_event <= cens),
      group,
      grid
    )
  }
  sims <- lapply(seq_len(n_sim), function(i) sim_once())
  # A failed replicate returns its error message instead of a data frame.
  failed <- vapply(sims, is.character, logical(1))
  if (any(failed)) {
    msg <- sprintf(
      "%d of %d replicates failed. First error, at %s",
      sum(failed),
      n_sim,
      sims[failed][[1]]
    )
    if (all(failed)) {
      stop(msg, call. = FALSE)
    }
    warning(msg, "\nThe failed replicates were dropped.", call. = FALSE)
  }
  sims <- do.call(rbind, sims[!failed])

  q <- stats::aggregate(
    surv ~ group + time,
    sims,
    function(s) stats::quantile(s, c(pi_probs[1], 0.5, pi_probs[2]))
  )
  sim <- data.frame(q[c("group", "time")], q$surv)
  names(sim)[3:5] <- c("lo", "median", "hi")

  # Censor the observed data at t_max too, so both curves end there.
  obs_status <- status * (time <= t_max)
  obs <- do.call(
    rbind,
    lapply(levels(group), function(g) {
      km <- survival::survfit(
        survival::Surv(cens, obs_status) ~ 1,
        subset = group == g
      )
      data.frame(
        group = g,
        time = c(0, km$time),
        surv = c(1, km$surv),
        censored = c(FALSE, km$n.censor > 0 & km$n.event == 0)
      )
    })
  )
  obs$group <- factor(obs$group, levels(group))
  sim$group <- factor(sim$group, levels(group))

  zph <- tryCatch(survival::cox.zph(fit), error = function(e) {
    warning(
      "cox.zph() failed, so the plot has no Schoenfeld p-value: ",
      conditionMessage(e),
      call. = FALSE
    )
    NULL
  })

  list(
    sim = sim,
    obs = obs,
    zph = zph,
    sample_size = sample_size,
    plot = .plot_vpc_cox(
      sim,
      obs,
      pi_probs,
      bootstrap,
      .zph_subtitle(zph, by),
      overlay
    )
  )
}

# TRUE when `data` (already without the rows the fit dropped) reproduces the
# fit's response and linear predictor row by row, i.e. it holds the fitted
# records in the fitted order.
.vpc_data_matches <- function(fit, data) {
  # Evaluate the Surv() response directly: model.frame(fit, data = data) looks
  # `data` up in the formula's environment and can find utils::data instead.
  tt <- stats::terms(fit)
  y_data <- tryCatch(
    eval(attr(tt, "variables")[[2L]], data, environment(tt)),
    error = function(e) NULL
  )
  if (is.null(y_data) || NROW(y_data) != NROW(fit$y)) return(FALSE)
  same_y <- isTRUE(all.equal(
    unclass(y_data), unclass(fit$y), check.attributes = FALSE
  ))
  lp <- tryCatch(
    stats::predict(fit, newdata = data, type = "lp"),
    error = function(e) NULL
  )
  same_y && !is.null(lp) &&
    isTRUE(all.equal(unname(lp), unname(fit$linear.predictors)))
}

.zph_subtitle <- function(zph, by) {
  if (is.null(zph)) {
    return(NULL)
  }
  p <- zph$table[, "p"]
  fmt_p <- function(x) format.pval(x, digits = 2, eps = 0.001)
  out <- sprintf(
    "Schoenfeld residual test of proportional hazards: global p = %s",
    fmt_p(p[["GLOBAL"]])
  )
  if (!is.null(by) && by %in% setdiff(names(p), "GLOBAL")) {
    out <- sprintf("%s; %s p = %s", out, by, fmt_p(p[[by]]))
  }
  out
}

# Returns a function that draws a censoring time for each subject who had the
# event, from the reverse Kaplan-Meier estimate G(t) = P(C > t) conditional on
# C exceeding their event time. Mass G leaves beyond the last time goes to
# t_max.
.reverse_km_sampler <- function(time, status, t_max) {
  rk <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  t_event <- time[status == 1]
  g_at_event <- c(1, rk$surv)[findInterval(t_event, rk$time) + 1]
  function() {
    target <- stats::runif(length(t_event)) * g_at_event
    # First time at which G falls to the target or below.
    j <- findInterval(-target, -rk$surv, left.open = TRUE) + 1
    out <- rep(t_max, length(t_event))
    out[j <= length(rk$time)] <- pmin(rk$time[j[j <= length(rk$time)]], t_max)
    out
  }
}

# Kaplan-Meier survival at each grid time, per group.
.km_at_grid <- function(time, event, group, grid) {
  do.call(
    rbind,
    lapply(levels(group), function(g) {
      km <- survival::survfit(
        survival::Surv(time, event) ~ 1,
        subset = group == g
      )
      s <- summary(km, times = grid, extend = TRUE)
      data.frame(group = g, time = s$time, surv = s$surv)
    })
  )
}

# Categorical colours for overlaid groups, used in this order (validated for
# colour-vision deficiency on a light background).
.vpc_group_colours <- c(
  "#2a78d6",
  "#eb6834",
  "#1baf7a",
  "#eda100",
  "#e87ba4",
  "#008300",
  "#6250d6",
  "#e34948"
)

.plot_vpc_cox <- function(
  sim,
  obs,
  pi_probs,
  bootstrap,
  subtitle = NULL,
  overlay = FALSE
) {
  ink <- "#1f1e1c"
  band <- "#3a6ea5"
  censored <- obs[obs$censored, ]
  p <- ggplot2::ggplot()
  if (overlay && nlevels(obs$group) > 1) {
    cols <- stats::setNames(
      .vpc_group_colours[seq_len(nlevels(obs$group))],
      levels(obs$group)
    )
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(time, ymin = lo, ymax = hi, fill = group),
        sim,
        alpha = 0.15
      ) +
      ggplot2::geom_line(
        ggplot2::aes(
          time,
          median,
          colour = group,
          linetype = "Simulated median"
        ),
        sim,
        linewidth = 0.6
      ) +
      ggplot2::geom_step(
        ggplot2::aes(time, surv, colour = group, linetype = "Observed"),
        obs,
        linewidth = 0.7
      ) +
      ggplot2::geom_point(
        ggplot2::aes(time, surv, colour = group),
        censored,
        shape = 3,
        size = 1.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        NULL,
        values = cols,
        aesthetics = c("colour", "fill")
      )
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(time, ymin = lo, ymax = hi),
        sim,
        fill = band,
        alpha = 0.25
      ) +
      ggplot2::geom_line(
        ggplot2::aes(time, median, linetype = "Simulated median"),
        sim,
        colour = band,
        linewidth = 0.6
      ) +
      ggplot2::geom_step(
        ggplot2::aes(time, surv, linetype = "Observed"),
        obs,
        colour = ink,
        linewidth = 0.6
      ) +
      ggplot2::geom_point(
        ggplot2::aes(time, surv),
        censored,
        shape = 3,
        size = 1.5,
        colour = ink
      )
    if (nlevels(obs$group) > 1) p <- p + ggplot2::facet_wrap(~group)
  }
  p +
    ggplot2::scale_linetype_manual(
      NULL,
      values = c("Observed" = "solid", "Simulated median" = "dashed")
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Time",
      y = "Proportion without event",
      subtitle = subtitle,
      caption = sprintf(
        "Shaded band: %g%%-%g%% prediction interval of simulated Kaplan-Meier curves%s; + censored",
        100 * pi_probs[1],
        100 * pi_probs[2],
        if (bootstrap) " (with parameter uncertainty, bootstrap)" else ""
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
}
