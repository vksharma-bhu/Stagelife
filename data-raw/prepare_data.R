## Builds the four real step-stress data sets bundled with the package.
## Run interactively with usethis::use_data(); not part of the installed
## package (data-raw/ is listed in .Rbuildignore).
##
## Naming convention: every bundled data set is named slt_<first author>_
## <year>, for consistency with every other exported object in this
## package (dslt/pslt/qslt/rslt, slt_families(), slt_fit objects, etc.
## all use the slt_/*slt* naming convention).
##
## Sources verified directly against the original papers (or, for
## slt_bobotas_kateri_2015 and slt_wang_fei_2003, against a later paper
## by this package's own author that correctly attributes them, which
## corrected two errors here: slt_kateri_nikolov_2024 (a real, but
## wrong, paper by a partially-overlapping author list) was replaced
## with the actual source, and slt_data3 -- previously shipped with no
## identified source at all -- was identified.

scale <- 100
slt_han_kundu_2014 <- list(
  y0 = c(0.140, 0.783, 1.324, 1.582, 1.716, 1.794, 1.883, 2.293, 2.660, 2.674,
         2.725, 3.085, 3.924, 4.396, 4.612, 4.892) * scale,
  y1 = c(5.002, 5.022, 5.082, 5.112, 5.147, 5.238, 5.244, 5.247, 5.305, 5.337,
         5.407, 5.408, 5.445, 5.483, 5.717) * scale,
  tau1 = 5 * scale,
  description = paste(
    "Simple step-stress data set (16 + 15 units) on the reliability of a",
    "solar lighting device, from Han & Kundu (2014), given here on a",
    "scale of 100x the original reported units (i.e. in hours, not",
    "hundreds of hours). n = 35 prototypes; temperature raised from 293K",
    "to 353K at the change-point; experiment terminated at 600 hours."
  ),
  source = paste(
    "Han, D. and Kundu, D. (2014). Inference for a step-stress model",
    "with competing risks for failure from the generalized exponential",
    "distribution under Type-I censoring. IEEE Transactions on",
    "Reliability, 64(1), 31-43."
  )
)

slt_bobotas_kateri_2015 <- list(
  y0 = c(8, 38, 72, 97, 122, 140, 163, 170, 188, 198, 223, 256, 257, 265, 448),
  y1 = c(608, 611, 614, 615, 616, 620, 623, 623, 624, 624, 631, 636, 646, 654,
         660, 673, 675, 680, 684, 692, 693, 730, 745),
  tau1 = 600,
  description = paste(
    "Simple step-stress data set (15 + 23 units) on the charge",
    "detrapping and dielectric breakdown of nanocrystalline embedded",
    "high-k devices, from Bobotas & Kateri (2015). Voltage (the stress",
    "factor) changed at 600 seconds; experiment terminated at 780 seconds."
  ),
  source = paste(
    "Bobotas, P. and Kateri, M. (2015). The step-stress tampered",
    "failure rate model under interval monitoring. Statistical",
    "Methodology, 27, 100-122."
  )
)

slt_wang_fei_2003 <- list(
  y0 = c(32, 54, 59, 86, 117, 123, 213, 267, 268, 273, 299, 311, 321, 333, 339,
         386, 408, 422, 435, 437, 476, 518, 570, 632, 666, 697, 796, 854, 858, 910),
  y1 = c(926, 929, 931, 946, 947, 973, 980, 985, 993, 1005, 1010, 1016, 1020,
         1023, 1026, 1045, 1046, 1059, 1082, 1096),
  tau1 = 910,
  description = paste(
    "Simple step-stress data set (30 + 20 units) on the reliability of",
    "electronic components, from Wang & Fei (2003). Operating",
    "temperature (the stress factor) raised from 100C to 150C at the",
    "change-point (910 hours); normal operating temperature is 25C."
  ),
  source = paste(
    "Wang, R. and Fei, H. (2003). Uniqueness of the maximum likelihood",
    "estimate of the Weibull distribution tampered failure rate model.",
    "Communications in Statistics - Theory and Methods, 32(12), 2321-2338."
  )
)

## Note: the source paper describes the nominal change-point as 96
## hours, but this is inconsistent with the recorded data itself: the
## smallest stage-1 (post-change) observation is 94.38 (Set-II), which
## would be impossible if the change-point were 96 (a unit cannot fail
## on stage 1 before the stage change happens). tau1 = 94 -- the only
## value consistent with both stages of the recorded data (all stage-0
## observations are <= 91.56; all stage-1 observations are >= 94.38) --
## is used here instead; see ?slt_zhu_2010.
slt_zhu_2010 <- list(
  y0 = c(12.07, 19.5, 22.1, 23.11, 24, 25.1, 26.9, 36.64, 44.1, 46.3, 54, 58.09,
         64.17, 72.25, 86.9, 90.09, 91.22,
         14, 17.95, 24, 26.46, 26.58, 28.06, 34, 36.13, 40.85, 41.11, 42.63,
         52.51, 62.68, 73.13, 83.63, 91.56),
  y1 = c(102.1, 105.1, 109.2, 114.4, 117.9, 121.9, 122.5, 123.6, 126.5, 130.1,
         94.38, 97.71, 101.53, 105.11, 112.11, 119.58, 120.2, 126.95, 129.25, 136.31),
  tau1 = 94,
  description = paste(
    "Simple step-stress data set (33 + 20 units, two sets of 32 light",
    "bulbs pooled) on bulb-filament fatigue, from Zhu (2010). Voltage",
    "(the stress factor) raised from 2.25V to 2.44V at the change-point;",
    "normal operating voltage is 2V. The source describes the nominal",
    "change-point as 96 hours, but the recorded data is only consistent",
    "with a change-point at or below 94.38 (the smallest stage-1",
    "observation); 94 is used here -- see Details in ?slt_zhu_2010."
  ),
  source = paste(
    "Zhu, Y. (2010). Optimal Design and Equivalency of Accelerated Life",
    "Testing Plans. PhD dissertation, Rutgers, The State University of",
    "New Jersey, School of Graduate Studies."
  )
)

usethis::use_data(slt_han_kundu_2014, overwrite = TRUE)
usethis::use_data(slt_bobotas_kateri_2015, overwrite = TRUE)
usethis::use_data(slt_wang_fei_2003, overwrite = TRUE)
usethis::use_data(slt_zhu_2010, overwrite = TRUE)
