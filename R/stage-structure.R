############################################################################
# Two-stage (step-stress / PALT) structural functions built from a single
# baseline distribution via the cumulative exposure model, plus the three
# identifiable parameter-sharing "cases" used for estimation.
############################################################################

#' Equivalent age under the elevated stress
#'
#' Under the cumulative exposure model (CEM), a unit that has accumulated
#' damage equivalent to surviving to time `tau` under stage-0 conditions
#' has, under stage-1 conditions, already "used up" an equivalent age
#' `nu`. This function computes `nu`.
#'
#' @details
#' The CEM's defining assumption is that only the *total accumulated
#' damage* determines a unit's remaining life, not the stress history
#' that produced it. Formalising "equal accumulated damage" as "equal
#' reliability", a unit surviving to \eqn{\tau} under the stage-0
#' distribution \eqn{F_0} is equivalent to a unit that has survived to
#' some age \eqn{\nu} under the stage-1 distribution \eqn{F_1}, where
#' \eqn{\nu} solves
#' \deqn{S_1(\nu;\theta_1) = S_0(\tau;\theta_0) \quad\Longleftrightarrow\quad
#'       F_1(\nu;\theta_1) = F_0(\tau;\theta_0),}
#' i.e.
#' \deqn{\nu = \nu(\tau;\theta_0,\theta_1) = Q_1\!\left(F_0(\tau;\theta_0);\ \theta_1\right),}
#' where \eqn{Q_1 = F_1^{-1}} is the stage-1 quantile function. This is
#' exactly `fam$quantile(fam$cdf(tau, theta0), theta1)` below: plug the
#' stage-0 cdf at `tau` into the stage-1 quantile function. \eqn{\nu}
#' is the building block for the population functions [pslt()]/[dslt()]/
#' [qslt()] and for the log-likelihood [loglik_slt()] -- see their
#' documentation for how it is used.
#'
#' @param family A family key (see [slt_families()]) or family object.
#' @param tau Non-negative scalar change-point time.
#' @param theta0 Numeric `c(lambda0, alpha0)`, baseline-stage parameters.
#' @param theta1 Numeric `c(lambda1, alpha1)`, elevated-stage parameters.
#' @return The equivalent age (a non-negative scalar).
#' @seealso [pslt()] for how `nu` combines with `tau1` into the
#'   population cdf; [loglik_slt()] for the resulting log-likelihood;
#'   [slt_tau()] for choosing `tau` as a stage-0 quantile.
#' @export
#' @examples
#' slt_equiv_age("weibull", tau = 5, theta0 = c(0.1, 1.5), theta1 = c(0.3, 1.5))
slt_equiv_age <- function(family, tau, theta0, theta1) {
  fam <- get_family(family)
  fam$quantile(fam$cdf(tau, theta0), theta1)
}

#' p-th quantile of the baseline-stage distribution
#'
#' Convenience wrapper typically used to set the change-point at a chosen
#' quantile of the stage-0 lifetime distribution, e.g.
#' `tau <- slt_tau(family, theta0, p = 0.5)` places the change-point at the
#' stage-0 median.
#'
#' @inheritParams slt_equiv_age
#' @param p Probability in `(0, 1)`.
#' @return The p-th quantile of the stage-0 distribution.
#' @export
#' @examples
#' slt_tau("weibull", theta0 = c(0.1, 1.5), p = 0.5)
slt_tau <- function(family, theta0, p = 0.5) {
  fam <- get_family(family)
  fam$quantile(p, theta0)
}

#' Population CDF/PDF/quantile function of a two-stage step-stress model
#'
#' `pslt()`, `dslt()` and `qslt()` implement the full stage life testing
#' design, in which every unit switches from the baseline to the
#' elevated stress at a common change-point `tau1`. Combined with the
#' cumulative-exposure equivalent age [slt_equiv_age()], the two
#' branches (pre-/post-tau1) are pieced together into a single,
#' continuous population CDF.
#'
#' @details
#' For a single unit switching from stage 0 to stage 1 at the fixed
#' time `tau1` (every unit switches -- the classical simple-step-stress
#' design, matching `prop_transfer = 1` in [rslt()]), the CEM
#' construction (see [slt_equiv_age()]) gives a population cdf
#' \deqn{F(t;\theta_0,\theta_1,\tau_1) = \begin{cases}
#'         F_0(t;\theta_0), & 0 < t \le \tau_1, \\
#'         F_1(t-\tau_1+\nu_1;\ \theta_1), & t > \tau_1,
#'       \end{cases}}
#' where \eqn{\nu_1 = \nu(\tau_1;\theta_0,\theta_1)} is the stage-1
#' equivalent age ([slt_equiv_age()]) at the moment of the switch.
#' Differentiating -- a constant \emph{shift} of the argument by
#' \eqn{\tau_1-\nu_1}, not a general change of variables, so no extra
#' Jacobian factor arises -- gives the matching population pdf
#' \deqn{f(t;\theta_0,\theta_1,\tau_1) = \begin{cases}
#'         f_0(t;\theta_0), & 0 < t \le \tau_1, \\
#'         f_1(t-\tau_1+\nu_1;\ \theta_1), & t > \tau_1.
#'       \end{cases}}
#' \eqn{F} is continuous at \eqn{\tau_1} by construction of \eqn{\nu_1}
#' (\eqn{F_1(\tau_1-\tau_1+\nu_1;\theta_1) = F_1(\nu_1;\theta_1) =
#' F_0(\tau_1;\theta_0)}, verified as an automatic consequence of the
#' equal-reliability definition of \eqn{\nu_1}, not a separate
#' assumption) -- but \eqn{F} is generally \emph{not} differentiable at
#' \eqn{\tau_1} (the hazard rate jumps there), an intrinsic, well-known
#' feature of the cumulative exposure model, not an artefact of this
#' implementation. `qslt()` inverts \eqn{F} directly (analytically, not
#' numerically): for \eqn{p \le F_0(\tau_1;\theta_0)} it returns
#' \eqn{F_0^{-1}(p;\theta_0)}; otherwise it returns \eqn{\tau_1-\nu_1+
#' F_1^{-1}(p;\theta_1)}, using each stage's own quantile function
#' (itself closed-form for every family in [slt_families()]).
#'
#' @inheritParams slt_equiv_age
#' @param t,x Numeric vector of times at which to evaluate the CDF/PDF.
#' @param p Numeric vector of probabilities at which to evaluate the
#'   quantile function.
#' @param tau1 Scalar change-point time (same argument as `tau` in
#'   [slt_equiv_age()]).
#' @return A numeric vector the same length as `t`/`x`/`p`.
#' @seealso [slt_equiv_age()] for `nu1`; [loglik_slt()] for the
#'   estimation log-likelihood; [rslt()] for simulating data under the
#'   actual random-withdrawal design; [fit_slt()] for maximum
#'   likelihood estimation.
#' @export
#' @examples
#' pslt(c(1, 6, 10), "weibull", c(0.1, 1.5), c(0.3, 1.5), tau1 = 5)
pslt <- function(t, family, theta0, theta1, tau1) {
  fam <- get_family(family)
  nu1 <- fam$quantile(fam$cdf(tau1, theta0), theta1)
  out <- numeric(length(t))
  out[is.na(t)] <- NA_real_
  ok <- !is.na(t)
  ind0 <- ok & t <= tau1
  ind1 <- ok & t > tau1
  out[ind0] <- fam$cdf(t[ind0], theta0)
  out[ind1] <- fam$cdf(t[ind1] - tau1 + nu1, theta1)
  out
}

#' @rdname pslt
#' @export
dslt <- function(x, family, theta0, theta1, tau1) {
  fam <- get_family(family)
  nu1 <- fam$quantile(fam$cdf(tau1, theta0), theta1)
  out <- numeric(length(x))
  out[is.na(x)] <- NA_real_
  ok <- !is.na(x)
  ind0 <- ok & x <= tau1
  ind1 <- ok & x > tau1
  out[ind0] <- fam$pdf(x[ind0], theta0)
  out[ind1] <- fam$pdf(x[ind1] - tau1 + nu1, theta1)
  out
}

#' @rdname pslt
#' @export
qslt <- function(p, family, theta0, theta1, tau1) {
  fam <- get_family(family)
  nu1 <- fam$quantile(fam$cdf(tau1, theta0), theta1)
  F_tau0 <- fam$cdf(tau1, theta0)
  out <- numeric(length(p))
  out[is.na(p)] <- NA_real_
  ok <- !is.na(p)
  ind0 <- ok & p <= F_tau0
  ind1 <- ok & p > F_tau0
  out[ind0] <- fam$quantile(p[ind0], theta0)
  out[ind1] <- tau1 - nu1 + fam$quantile(p[ind1], theta1)
  out
}

############################################################################
## Parameter-sharing "cases"
##
## Four logically possible sharing structures exist for the pair
## theta0 = (lambda0, alpha0), theta1 = (lambda1, alpha1):
##
##   full     : (lambda0, alpha0), (lambda1, alpha1)   -- 4 free parameters
##   common   : (lambda0, alpha ), (lambda1, alpha )   -- 3 free parameters
##   shape1   : (lambda0, alpha0), (lambda1, 1     )   -- 3 free parameters
##   shape11  : (lambda0, 1     ), (lambda1, 1     )   -- 2 free parameters
##
## 'full' is NOT offered by this package: with a single change-point and
## no replicated stress levels, the likelihood surface in (alpha0, alpha1)
## is essentially flat along a ridge (independent shape parameters at each
## stage are not separably identifiable from a two-stage design), which in
## simulation manifests as a high rate of non-convergence / singular
## Hessians and estimates that are extremely sensitive to starting values.
## 'shape1' fixes alpha1 = 1 and is the structure used for the primary
## analysis in the source study; 'shape11' additionally fixes alpha0 = 1
## and is the natural nested null used for the likelihood-ratio "shape = 1"
## test; 'common' (a shared shape across stages, both scales free) is
## offered as a useful, well-identified alternative not considered in the
## original study.
############################################################################

.slt_cases <- list(
  common = list(
    key = "common", npar = 3L,
    par_names = c("lambda0", "lambda1", "alpha"),
    label = "theta0=(lambda0,alpha), theta1=(lambda1,alpha)  [shared shape]"
  ),
  shape1 = list(
    key = "shape1", npar = 3L,
    par_names = c("lambda0", "alpha0", "lambda1"),
    label = "theta0=(lambda0,alpha0), theta1=(lambda1,1)  [stage-1 shape fixed at 1]"
  ),
  shape11 = list(
    key = "shape11", npar = 2L,
    par_names = c("lambda0", "lambda1"),
    label = "theta0=(lambda0,1), theta1=(lambda1,1)  [both shapes fixed at 1]"
  )
)

#' Identifiable parameter-sharing cases for two-stage estimation
#'
#' Describes the three parameter-sharing structures between the
#' baseline-stage parameters `theta0 = (lambda0, alpha0)` and the
#' elevated-stage parameters `theta1 = (lambda1, alpha1)` supported for
#' estimation. A fourth, fully unrestricted structure (`alpha0` and
#' `alpha1` both free) is deliberately not supported; see Details.
#'
#' @details
#' Four parameter-sharing structures are logically possible for the pair
#' \eqn{\theta_0=(\lambda_0,\alpha_0)}, \eqn{\theta_1=(\lambda_1,\alpha_1)}:
#' \deqn{\begin{array}{lll}
#'   \texttt{full}    : & \theta_0=(\lambda_0,\alpha_0), & \theta_1=(\lambda_1,\alpha_1)
#'                        \quad\text{(4 free parameters)} \\
#'   \texttt{common}  : & \theta_0=(\lambda_0,\alpha),   & \theta_1=(\lambda_1,\alpha)
#'                        \quad\text{(3 free parameters, shared shape)} \\
#'   \texttt{shape1}  : & \theta_0=(\lambda_0,\alpha_0), & \theta_1=(\lambda_1,1)
#'                        \quad\text{(3 free parameters)} \\
#'   \texttt{shape11} : & \theta_0=(\lambda_0,1),        & \theta_1=(\lambda_1,1)
#'                        \quad\text{(2 free parameters)}
#' \end{array}}
#' This package supports the last three; \code{full} is deliberately
#' **not** offered (see "Why `full` is excluded" below).
#'
#' `shape1`/`shape11` fix the elevated-stage shape to the *numeric
#' value* 1 -- this is a structural constraint, not a claim that the
#' elevated stage becomes Exponential. The two only coincide for the
#' four families whose `alpha = 1` boundary happens to be Exponential
#' (`weibull`, `expexponential`, `moexponential`, `gamma`; see
#' `?slt_dist`); for the other eight families, `alpha1 = 1` is simply
#' that family's own (otherwise unnamed) member at shape 1 -- e.g. Chen
#' with its shape fixed at 1, not an Exponential elevated stage.
#' Fixing the shape is still meaningful regardless: it
#' reduces the free parameter count the same way for every family (see
#' "Why `full` is excluded" below), and "is the elevated-stage shape
#' exactly 1?" is a concrete, testable nested hypothesis via
#' [slt_lrt()] whether or not that constraint has a textbook name.
#'
#' Each case defines a bijection between its free-parameter vector `par`
#' (on the natural, positive scale, in the order given by
#' `case_info(case)$par_names`) and `(theta0, theta1)`, implemented by
#' [slt_expand_params()] and inverted by [slt_reduce_params()]:
#' \deqn{\texttt{common}:\ (\lambda_0,\lambda_1,\alpha) \mapsto
#'         \big((\lambda_0,\alpha),\ (\lambda_1,\alpha)\big),}
#' \deqn{\texttt{shape1}:\ (\lambda_0,\alpha_0,\lambda_1) \mapsto
#'         \big((\lambda_0,\alpha_0),\ (\lambda_1,1)\big),}
#' \deqn{\texttt{shape11}:\ (\lambda_0,\lambda_1) \mapsto
#'         \big((\lambda_0,1),\ (\lambda_1,1)\big).}
#' `fit_slt()` optimises the log-likelihood ([loglik_slt()], via
#' [neg_loglik_slt()]) over `log(par)` in this order, then reports the
#' fitted `par` under exactly these names -- so, for instance, `case =
#' "shape1"`'s second free parameter is always labelled `alpha0` and is
#' always the baseline-stage shape, never accidentally the elevated-stage
#' rate or any other quantity. This parameter-index bookkeeping (`
#' case_info()$par_names`, [slt_expand_params()], [slt_reduce_params()],
#' and `fit_slt()`'s internal default-initial-value heuristic all sharing
#' one consistent order) is exactly the kind of place a silent
#' mislabelling bug could hide, so it has been checked directly: the
#' round trip `slt_expand_params(slt_reduce_params(theta0, theta1,
#' case), case)` reproduces `(theta0, theta1)` exactly for all three
#' cases, and a full simulate-then-refit check (large `n`, known
#' `theta0`/`theta1`) recovers the true values under the correct
#' parameter names for all three cases (see the package's test suite,
#' `test-fit.R`).
#'
#' `full` (both shapes free) is excluded because, with a single
#' change-point and no replicated stress levels, its extra free
#' parameter is empirically non-identifiable (a flat likelihood ridge,
#' frequent non-convergence); `shape1`/`shape11` remove this by fixing
#' one or both shapes to 1, `common` by tying them together instead.
#' `shape11` is nested in both `common` and `shape1` (the natural null
#' model for [slt_lrt()]); `shape1` and `common` are not nested in each
#' other, so compare those two with [slt_model_comparison()] instead.
#'
#' @return A data.frame with columns `case`, `n_par`, `parameters`,
#'   `description`.
#' @seealso [slt_expand_params()] and [slt_reduce_params()] for the
#'   parameter mapping implementing each case; [fit_slt()] for
#'   estimation; [slt_lrt()] for the nested likelihood-ratio test;
#'   [slt_model_comparison()] for comparing non-nested cases.
#' @export
#' @examples
#' slt_cases()
slt_cases <- function() {
  data.frame(
    case = vapply(.slt_cases, `[[`, character(1), "key"),
    n_par = vapply(.slt_cases, `[[`, integer(1), "npar"),
    parameters = vapply(.slt_cases, function(z) paste(z$par_names, collapse = ", "), character(1)),
    description = vapply(.slt_cases, `[[`, character(1), "label"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
case_info <- function(case) {
  case <- match.arg(case, names(.slt_cases))
  .slt_cases[[case]]
}

#' Expand a case-specific parameter vector into `(theta0, theta1)`
#'
#' @param par Numeric vector on the **natural** (positive) scale, with
#'   length and order given by `case_info(case)$par_names`.
#' @param case One of `"common"`, `"shape1"`, `"shape11"`.
#' @return A list with numeric elements `theta0` and `theta1`, each of
#'   length 2.
#' @keywords internal
slt_expand_params <- function(par, case = c("common", "shape1", "shape11")) {
  case <- match.arg(case)
  switch(case,
    common  = list(theta0 = c(par[1], par[3]), theta1 = c(par[2], par[3])),
    shape1  = list(theta0 = c(par[1], par[2]), theta1 = c(par[3], 1)),
    shape11 = list(theta0 = c(par[1], 1),      theta1 = c(par[2], 1))
  )
}

#' Reduce full `(theta0, theta1)` parameters to a case's free parameters
#'
#' The inverse of [slt_expand_params()]; mostly useful for extracting the
#' "true" parameter vector matching a given case in simulation studies.
#' For `"common"`, `theta0[2]` (rather than an average) is used as the
#' shared shape and a warning is issued if `theta0[2] != theta1[2]`.
#'
#' @param theta0,theta1 Numeric `c(lambda, alpha)` vectors.
#' @param case One of `"common"`, `"shape1"`, `"shape11"`.
#' @return A named numeric vector matching `case_info(case)$par_names`.
#' @export
#' @examples
#' slt_reduce_params(c(0.2, 1.5), c(0.4, 1.5), case = "common")
slt_reduce_params <- function(theta0, theta1, case = c("common", "shape1", "shape11")) {
  case <- match.arg(case)
  out <- switch(case,
    common = {
      if (isTRUE(all.equal(theta0[2], theta1[2]))) {
        c(lambda0 = theta0[1], lambda1 = theta1[1], alpha = theta0[2])
      } else {
        warning("theta0[2] != theta1[2]; using theta0[2] as the shared shape.",
                call. = FALSE)
        c(lambda0 = theta0[1], lambda1 = theta1[1], alpha = theta0[2])
      }
    },
    shape1 = c(lambda0 = theta0[1], alpha0 = theta0[2], lambda1 = theta1[1]),
    shape11 = c(lambda0 = theta0[1], lambda1 = theta1[1])
  )
  out
}
