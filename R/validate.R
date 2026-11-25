# Argument validators shared by the exported functions. Each stops with a
# message naming the argument, e.g. "`n_boot` must be ...".

.check_count <- function(x, name, min = 1L) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok)
    stop(sprintf("`%s` must be a single whole number of at least %d.",
                 name, min), call. = FALSE)
  as.integer(x)
}

.check_prob <- function(x, name) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
  if (!ok)
    stop(sprintf("`%s` must be a single number strictly between 0 and 1.",
                 name), call. = FALSE)
  invisible(x)
}

.check_seed <- function(x, name = "seed") {
  if (is.null(x)) return(invisible(x))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(sprintf("`%s` must be a single number or NULL.", name), call. = FALSE)
  invisible(x)
}

# Numeric vector of finite values in [lower, upper]; NULL allowed.
.check_range <- function(x, name, lower = -Inf, upper = Inf) {
  if (is.null(x)) return(invisible(x))
  ok <- is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
    all(x >= lower) && all(x <= upper)
  if (!ok)
    stop(sprintf("`%s` must be finite numbers between %s and %s.",
                 name, format(lower), format(upper)), call. = FALSE)
  invisible(x)
}
