############################################################################
# Four real simple step-stress data sets, each a list with elements
# `y0` (stage-0 failure times), `y1` (stage-1 failure times), `tau1`
# (the change-point), `description`, and `source`. Naming convention:
# slt_<first author>_<year>, matching every other exported object in
# this package.
############################################################################

#' Han & Kundu (2014) simple step-stress data
#'
#' A simple step-stress life-testing data set on the reliability of a
#' solar lighting device, with 16 units observed under the baseline
#' stress (normal operating temperature, 293K) and 15 under the elevated
#' stress (353K) after a common change-point at `tau1 = 500` (hours).
#'
#' @details
#' The original experiment placed `n = 35` units on test, but only
#' these 31 failure times are reported in the source paper -- the
#' remaining 4 units did not fail before the experiment ended. This
#' package's model has no right-censoring mechanism (see
#' `?stagelife-package`), so those 4 units are not, and cannot be,
#' represented here.
#'
#' @format A list with elements:
#' \describe{
#'   \item{y0}{Numeric vector of 16 stage-0 (baseline stress) failure times.}
#'   \item{y1}{Numeric vector of 15 stage-1 (elevated stress) failure times.}
#'   \item{tau1}{The change-point time.}
#'   \item{description, source}{Character strings.}
#' }
#' @source Han, D. and Kundu, D. (2014). Inference for a step-stress
#'   model with competing risks for failure from the generalized
#'   exponential distribution under Type-I censoring. *IEEE Transactions
#'   on Reliability*, 64(1), 31-43.
#' @seealso [slt_bobotas_kateri_2015], [slt_wang_fei_2003],
#'   [slt_zhu_2010] for the other three bundled data sets; [fit_slt()]
#'   for fitting a model to data in this format.
#' @examples
#' data(slt_han_kundu_2014)
#' fit <- fit_slt(list(Y = slt_han_kundu_2014$y0, Z = slt_han_kundu_2014$y1),
#'                   "weibull", "shape1", tau1 = slt_han_kundu_2014$tau1)
#' summary(fit)
"slt_han_kundu_2014"

#' Bobotas & Kateri (2015) simple step-stress data
#'
#' A simple step-stress life-testing data set on the charge detrapping
#' and dielectric breakdown of nanocrystalline embedded high-k devices,
#' with 15 units observed under the baseline stress and 23 under the
#' elevated stress after a common change-point at `tau1 = 600` (seconds).
#'
#' @format A list with elements `y0`, `y1`, `tau1`, `description`, `source`
#'   (see [slt_han_kundu_2014] for the structure).
#' @source Bobotas, P. and Kateri, M. (2015). The step-stress tampered
#'   failure rate model under interval monitoring. *Statistical
#'   Methodology*, 27, 100-122.
#' @seealso [slt_han_kundu_2014], [slt_wang_fei_2003], [slt_zhu_2010].
#' @examples
#' data(slt_bobotas_kateri_2015)
#' str(slt_bobotas_kateri_2015)
"slt_bobotas_kateri_2015"

#' Wang & Fei (2003) simple step-stress data
#'
#' A simple step-stress life-testing data set on the reliability of
#' electronic components under temperature stress, with 30 units
#' observed under the baseline stress and 20 under the elevated stress
#' after a common change-point at `tau1 = 910` (hours).
#'
#' @format A list with elements `y0`, `y1`, `tau1`, `description`, `source`
#'   (see [slt_han_kundu_2014] for the structure).
#' @source Wang, R. and Fei, H. (2003). Uniqueness of the maximum
#'   likelihood estimate of the Weibull distribution tampered failure
#'   rate model. *Communications in Statistics - Theory and Methods*,
#'   32(12), 2321-2338.
#' @seealso [slt_han_kundu_2014], [slt_bobotas_kateri_2015], [slt_zhu_2010].
#' @examples
#' data(slt_wang_fei_2003)
#' str(slt_wang_fei_2003)
"slt_wang_fei_2003"

#' Zhu (2010) simple step-stress data
#'
#' A simple step-stress life-testing data set on light-bulb-filament
#' fatigue (two sets of 32 bulbs pooled), with 33 units observed under
#' the baseline stress and 20 under the elevated stress after a common
#' change-point at `tau1 = 94`.
#'
#' @details
#' The original experiment used two sets of 32 bulbs (`n = 64` total),
#' but only these 53 failure times are bundled here -- the remaining 11
#' units did not fail before their set's test ended. This package's
#' model has no right-censoring mechanism (see `?stagelife-package`),
#' so those 11 units are not, and cannot be, represented here.
#'
#' The source describes the nominal change-point as 96 hours, but this
#' is inconsistent with the recorded data itself: the smallest stage-1
#' (post-change) observation is 94.38 (Set-II), which would be
#' impossible if the change-point were 96 (a unit cannot fail on stage 1
#' before the stage change happens). `tau1 = 94` -- the only value
#' consistent with both stages of the recorded data (every stage-0
#' observation is `<= 91.56`; every stage-1 observation is `>= 94.38`)
#' -- is used here instead.
#'
#' @format A list with elements `y0`, `y1`, `tau1`, `description`, `source`
#'   (see [slt_han_kundu_2014] for the structure).
#' @source Zhu, Y. (2010). *Optimal Design and Equivalency of
#'   Accelerated Life Testing Plans*. PhD dissertation, Rutgers, The
#'   State University of New Jersey, School of Graduate Studies.
#' @seealso [slt_han_kundu_2014], [slt_bobotas_kateri_2015], [slt_wang_fei_2003].
#' @examples
#' data(slt_zhu_2010)
#' str(slt_zhu_2010)
"slt_zhu_2010"
