############################################################################
# Data generation: the stage life testing (SLT) model of Laumen & Cramer
# (2019, Procedure 2.3) is this package's only data-generating process.
# At the stage-change time tau1, R*_1 of the units still surviving are
# moved to the elevated stage (a fixed target proportion, "typeP", or a
# fixed target count, "typeM"); the rest continue failing under the
# original stage conditions. The classical "simple step-stress" design (every
# survivor switches) is simply the pi1 = 1 special case of "typeP" --
# Laumen & Cramer's Remark 2.2(4) -- so no separate function is needed
# for it.
############################################################################

#' Simulate data from the two-stage life-testing model
#'
#' Simulates data from the two-stage life-testing model of Laumen &
#' Cramer (2019, Procedure 2.3).
#'
#' @details
#' At the stage-change time `tau1`, `R*_1` of the units that have *not
#' yet failed* are randomly selected and moved to the elevated stage
#' (cumulative-exposure shifted); the remaining survivors continue
#' failing under the original stage-0 conditions and are never observed
#' under stage 1.
#'
#' Two rules for generating `R*_1` are supported, following the paper:
#' \describe{
#'   \item{`"typeP"`}{A fixed target *proportion* `pi1` of the survivors
#'     is transferred: `R*_1 = floor(pi1 * (n - D1))`, where `D1` is the
#'     (random) number of stage-0 failures observed before `tau1`.
#'     `pi1 = 1` transfers *every* survivor, which is exactly the
#'     classical "simple step-stress" design (Remark 2.2(4) of the
#'     paper) -- no separate function is needed for that special case.}
#'   \item{`"typeM"`}{A fixed target *count* `R0_1` is transferred, capped
#'     by the number of available survivors: `R*_1 = min(n - D1, R0_1)`.}
#' }
#'
#' @inheritParams slt_equiv_age
#' @param n Number of units placed on test.
#' @param tau1 Time at which the transfer decision is made.
#' @param design `"typeP"` (default, fixed target proportion) or
#'   `"typeM"` (fixed target count); see Details.
#' @param prop_transfer Proportion (in `[0, 1]`) of survivors at `tau1`
#'   transferred to the elevated stage, for `design = "typeP"`. Set to
#'   `1` for the classical simple-step-stress design.
#' @param R0_1 Target number of survivors transferred (capped by the
#'   number actually available), for `design = "typeM"`.
#' @return A list with `Y` (all stage-0 failure times: early failures
#'   before `tau1` plus the non-transferred survivors, who keep failing
#'   under stage-0), `Z` (stage-1 observed failure times for the
#'   transferred units), `n_early` (`D1`, number failing before `tau1`)
#'   and `n_transferred` (`R*_1`).
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(100, "weibull", theta0 = c(0.05, 1.5), theta1 = c(0.2, 1.5),
#'           tau1 = 5, design = "typeP", prop_transfer = 0.5)
#' str(d)
#'
#' d_m <- rslt(100, "weibull", theta0 = c(0.05, 1.5), theta1 = c(0.2, 1.5),
#'             tau1 = 5, design = "typeM", R0_1 = 30)
#' str(d_m)
#'
#' # The classical simple-step-stress design (every survivor transferred):
#' d_full <- rslt(100, "weibull", theta0 = c(0.05, 1.5), theta1 = c(0.2, 1.5),
#'                tau1 = 5, design = "typeP", prop_transfer = 1)
#' str(d_full)
rslt <- function(n, family, theta0, theta1, tau1,
                  design = c("typeP", "typeM"),
                  prop_transfer = 0.5, R0_1 = NULL) {
  design <- match.arg(design)
  check_positive(theta0, "theta0"); check_positive(theta1, "theta1")
  if (tau1 <= 0) stop("tau1 must be > 0.", call. = FALSE)
  if (design == "typeP" && (prop_transfer < 0 || prop_transfer > 1)) {
    stop("prop_transfer must be in [0, 1] for design = 'typeP'.", call. = FALSE)
  }
  if (design == "typeM" && (is.null(R0_1) || R0_1 < 0 || R0_1 != round(R0_1))) {
    stop("R0_1 must be a non-negative integer for design = 'typeM'.", call. = FALSE)
  }
  fam <- get_family(family)
  U <- sort(stats::runif(n))
  V0 <- fam$quantile(U, theta0)
  nu1 <- slt_equiv_age(family, tau1, theta0, theta1)

  early <- which(V0 <= tau1)
  survivors <- setdiff(seq_len(n), early)
  n_transfer <- if (design == "typeP") {
    floor(prop_transfer * length(survivors))
  } else {
    min(length(survivors), R0_1)
  }
  transferred <- if (n_transfer > 0L && n_transfer <= length(survivors)) {
    sort(sample(survivors, n_transfer))
  } else integer(0L)
  stayed <- setdiff(survivors, transferred)

  V1_transferred <- if (length(transferred)) {
    fam$quantile(U[transferred], theta1) - nu1 + tau1
  } else numeric(0L)

  list(Y = sort(c(V0[early], V0[stayed])),
       Z = sort(V1_transferred),
       n_early = length(early),
       n_transferred = length(transferred))
}
