#' @keywords internal
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

#' Guarded matrix inverse for variance-covariance matrices
#'
#' Attempts to invert a Hessian and returns `NULL` (rather than throwing an
#' error) if the matrix is singular, not positive-definite, or otherwise
#' unusable as a covariance matrix.
#'
#' @param hessian A numeric, symmetric matrix (typically the Hessian of the
#'   negative log-likelihood at the MLE).
#' @return The inverse matrix, or `NULL` if it is not a valid covariance
#'   matrix (any non-finite or non-positive diagonal element).
#' @keywords internal
safe_vcov <- function(hessian) {
  if (is.null(hessian) || any(!is.finite(hessian))) return(NULL)
  vc <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(vc)) return(NULL)
  d <- diag(vc)
  if (any(!is.finite(d)) || any(d <= 0)) return(NULL)
  # symmetrize to guard against numerical asymmetry
  (vc + t(vc)) / 2
}

#' Validate a numeric parameter vector is strictly positive and finite
#' @keywords internal
check_positive <- function(x, name = "theta") {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf("All elements of '%s' must be finite and strictly positive.", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a confidence level is a single number strictly between 0 and 1
#'
#' Guards against a silent, structurally-invalid confidence interval
#' (lower bound > upper bound) being returned with no warning for
#' `level` outside `(0, 1)` -- e.g. a fat-fingered `level = 95` instead
#' of `level = 0.95`. Found during an audit to be missing (with exactly
#' this silent-failure consequence) from every CI-producing function in
#' this package; now called by all of them.
#' @keywords internal
check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1 (got ",
         format(level), ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Round a numeric vector/data.frame column set, silently ignoring non-numeric
#' @keywords internal
round_df <- function(df, digits = 4) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  df
}

#' Identify a failed Monte Carlo / bootstrap replicate
#'
#' A replicate is considered failed if its worker function returned `NULL`
#' (the convention used throughout this package for "this replicate could
#' not be computed"), or if it is an error/condition object. The latter
#' case arises specifically when [parallel::mclapply()] is used: a forked
#' worker that hits an uncaught error does not propagate the error to the
#' caller, it returns a `"try-error"` object in that slot instead, which
#' would otherwise silently corrupt downstream aggregation (e.g.
#' `rbind()`-ing a character error message in among numeric estimates).
#' Every parallel replicate loop in this package (bootstrap, goodness of
#' fit, simulation studies) filters its results through
#' [drop_failed_replicates()] for this reason, on top of wrapping each
#' individual replicate in its own `tryCatch()` so a single numerically
#' fragile draw (most commonly seen with the Power Lindley family; see
#' `vignette("stagelife-intro")`) is skipped rather than aborting the
#' whole run.
#'
#' @param x A single replicate result.
#' @return A logical scalar.
#' @keywords internal
is_failed_replicate <- function(x) {
  is.null(x) || inherits(x, "condition") || inherits(x, "try-error")
}

#' @rdname is_failed_replicate
#' @param res A list of replicate results (e.g. from [lapply()] or
#'   [parallel::mclapply()]).
#' @return `drop_failed_replicates()` returns the subset of `res` that are
#'   not failures.
#' @keywords internal
drop_failed_replicates <- function(res) {
  Filter(Negate(is_failed_replicate), res)
}
