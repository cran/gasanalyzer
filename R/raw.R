#' Get user and factory calibration factors for use in calculations.
#'
#' This reshapes the calibration factors stored in list columns (passed
#' as arguments) in a matrix that can more easily be used in equations.
#' Only calibration factors used in calculation are kept and the list elements
#' are returned as columns in the output matrix.
#'
#' @param headcal a list-column containing the instrument head calibration
#' @param factcal a list-column containing the instrument factory calibration
#'
#' @returns a character matrix with calibration factors (in cols) for every
#'   rows
#'
#' @noRd
get_cal <- function(headcal, factcal) {

  # things we need
  need_fc <- get("needed_cals", .gasanalyzerEnv)[["fc"]]
  need_hc <-  get("needed_cals", .gasanalyzerEnv)[["hc"]]

  #prevent bugs:
  if (length(headcal) == 0L) headcal <- list(co2aspan1 = NA)
  if (length(factcal) == 0L) headcal <- list(caldata = NA)
  out <- matrix(NA_character_, nrow = length(need_hc) + length(need_fc),
                ncol = length(headcal),
                dimnames = list(c(need_hc, need_fc)))

  hcmat <- unlist(headcal)
  hcmat <- hcmat[names(hcmat) %in% need_hc]

  if (anyNA(hcmat))
      warning("Missing head calibration for some data!")

  out[need_hc, ] <- hcmat[need_hc]

  fcmat <- unlist(factcal)
  fcmat <- fcmat[names(fcmat) %in% need_fc]

  if (anyNA(factcal) && getOption("gasanalyzer.calibration.warning", TRUE)) {
    warning("Missing factory calibration data. This will lead ",
            "to NA values if gas mol fractions are recalculated. ",
            "No further warnings will be displayed!")
    options("gasanalyzer.calibration.warning" = FALSE)
  }

  #this if shouldnt be needed:
  if (length(out[need_fc, ]) / ncol(out) == length(fcmat[need_fc]))
    out[need_fc, ] <- fcmat[need_fc]
  else
    out[need_fc, ] <- NA

  t(out)
}

#' Calculate the oxygen correction factors for water
#'
#' This uses the O2 concentration in ``%` to calculate the correction factor
#' for the water mol fraction using the "bb" (band-broadening) coefficients
#' stored in the factory calibration data provided by headcal
#'
#' @param hc a matrix with the same number of rows as O2, providing
#'   bb.hx0 and bb.hx1 cal factors in one of the rows.
#' @param O2 a vector with the oxygen concentration in `%`.
#' @returns a correction factor for the water mol fraction
#'
#' @noRd
psiH2O <- function(hc, O2) {
  as.numeric(hc[ ,"bb.hx0"]) + as.numeric(hc[ ,"bb.hx1"]) * O2@"%"@.
}

#' Calculate the oxygen correction factors for CO2
#'
#' This uses the O2 concentration in `%` and water mol fraction in mmol / mol
#' to calculate the correction factor for the CO2 using the bb
#' (band-broadening) coefficients stored in the factory calibration data
#' provided by headcal
#'
#' @param hc a matrix with the same number of rows as O2, providing
#'   bb.ch and bb.cx cal factors in one of the rows.
#' @param O2 a vector with the oxygen concentration in `%`.
#' @param H2O a vector with the water mol fraction in mmol / mol.
#'
#' @return a correction factor for the water mol fraction
#' @noRd
psiCO2 <- function(hc, O2, H2O) {
  1 + (as.numeric(hc[, "bb.ch"]) - 1) * H2O@"mmol*mol^-1"@. / 1000 +
    (as.numeric(hc[, "bb.cx"]) - 1) * O2@"%"@. / 100
}

#' Apply span and pressure corrections to absorptance raw values
#'
#' This is a convenience function to calculate span corrected absorptance
#' values from raw absorptance values for the li6800. It corrects for span,
#' span drift and pressure.
#'
#' @param hc a matrix with the same number of rows as ab, providing
#'   span1, span2 and pcorr_ for the relevant h2o and co2 channels in one of the
#'   rows.
#' @param gas a character string with gas, either co2 or h2o
#' @param ch a character string with the channel, either a (sample) or b (ref)
#' @param pa a vector with ambient pressure (Meas.Pa)
#' @param ab a vector with raw absorptance values (Raw.abs_c_a) matching the
#'   gas and channel specified in the previous arguments
#'
#' @return a vector with span and pressure corrected absorptance values
#' @noRd
span_abs <-function(hc, gas = c("co2", "h2o"), ch = c("a", "b"), pa, ab) {
  pa <- pa@"kPa"@.
  gas <- match.arg(gas, choices = c("co2", "h2o"))
  ch <- match.arg(ch, choices = c("a", "b"))
  fc <- paste0("irga_", ch, ".", gas, ".", c("s1", "s2", "s3"))
  #FIXME: compat for new methods
  if(any(grepl("hc.",colnames(hc)))) cmpt <- "hc." else cmpt <-""
  span1 <- paste0(cmpt, gas, ch, "span1")
  span2 <- paste0(cmpt, gas, ch, "span2")
  pcor <- paste0("bb.pcorr_", substr(gas, 1, 1), c(1,2,3))
  x <- 1 - pa / 98
  #span drift:
  #ab <- ab * (1 + Tirga - as.numeric(hc[, fc[3]])) *
  #    (as.numeric(hc[, fc[1]]) + as.numeric(hc[, fc[2]]) * ab)

  (ab * (as.numeric(hc[, span1]) + ab * as.numeric(hc[, span2])) /
      (pa * (1 + as.numeric(hc[, pcor[1]]) * x +
               (as.numeric(hc[, pcor[2]]) + as.numeric(hc[, pcor[3]]) * ab)
             * x * x)))@"kPa^-1"
}

#' Calculate mol fraction/K from absorptance/P
#'
#' This is a convenience function to evaluate the calibration polynomial
#' defined in the factory calibration for the li6800. It relates the span
#' corrected absorptance per unit pressure (adjusted with psiCO2 for CO2)
#' to a mol fraction per unit temperature. The returned value may need further
#' adjustments for O2 and H2O (band-broadening).
#'
#' @param hc a matrix with the same number of rows as x, providing
#'   span1 and span2 for the relevant h2o and co2 channel in one of the rows.
#' @param gas a character string with gas, either co2 or h2o
#' @param ch a character string with the channel, either a (sample) or b (ref)
#' @param x a vector with span corrected absorptance values per unit pressure
#'  (the output of span divided through pressure, or psiCO2 * pressure)
#'
#' @return a vector with mol fraction per K for the indicated gas and channel
#'
#' @noRd
abs2frac <-function(hc, gas = c("co2", "h2o"), ch = c("a", "b"), x) {
  gas <- match.arg(gas, choices = c("co2", "h2o"))
  ch <- match.arg(ch, choices = c("a", "b"))
  fc <- paste0("irga_", ch, ".", gas, ".", c("a1", "a2", "a3", "a4", "a5"))

  x <- x@"kPa^-1"@.
  if (gas == "co2") {
    hc <- hc[, fc]
    #as.numeric would change to vector:
    storage.mode(hc) <- "numeric"
    (hc[, 1] * x + hc[, 2] * x^2 + hc[, 3] * x^3 + hc[, 4] * x^4 +
        hc[, 5] * x^5)@"\U00B5mol*mol^-1*K^-1"
  }
  else {
    hc <- hc[, fc[1:3]]
    #as.numeric would change to vector:
    storage.mode(hc) <- "numeric"
    (hc[, 1] * x + hc[, 2] * x^2 + hc[, 3] * x^3)@"mmol*mol^-1*K^-1"
  }
}

#' Calculate absorptance/P from mol fraction/K
#'
#' This is a convenience function to back calculate absorptance/P from
#' mol fractions/L using the calibration polynomials defined in the factory
#' calibration for the li6800. For CO2, the returned value is convenient for
#' applying O2 corrections with psiCO2, and then recalculating mol fraction/K
#' using abs2frac.
#'
#' @param hc a matrix with the same number of rows as fr, providing
#'   span1 and span2 for the relevant h2o and co2 channel in one of the rows.
#' @param gas a character string with gas, either co2 or h2o
#' @param ch a character string with the channel, either a (sample) or b (ref)
#' @param fr a vector with a mol fraction per K for the indicated gas and
#'   channel
#'
#' @return a vector a span corrected absorptance per unit pressure
#' @noRd
frac2abs <- function(hc, gas = c("co2", "h2o"), ch = c("a", "b"), fr) {
  gas <- match.arg(gas, choices = c("co2", "h2o"))
  ch <- match.arg(ch, choices = c("a", "b"))
  fc <- paste0("irga_", ch, ".", gas, ".", c("a1", "a2", "a3", "a4", "a5"))
  # I've some R glue for a cpp header that does real root isolation faster,
  # but polyroot performs well for low order polynomials so I use that for now

  if (gas == "co2") {
    apply(cbind(-fr@"umol/mol/K"@., matrix(as.numeric(hc[, fc]),
                                           nrow=length(fr))), 1,
          function(a) { if(any(is.na(a))) NA else {
            y <- polyroot(a);
            min(Re(y)[abs(Im(y)) < 1e-12])
          }
          })@"kPa^-1"
  } else {
    apply(cbind(-fr@"mmol/mol/K"@., matrix(as.numeric(hc[, fc[1:3]]),
                                           nrow=length(fr))), 1,
          function(a) { if(any(is.na(a))) NA else {
            y <- polyroot(a);
            min(Re(y)[abs(Im(y)) < 1e-12])
          }
          })@"kPa^-1"
  }

}

#' Calculate CO2 and H2O values that are uncorrected for O2
#'
#' The 6800 stores the O2 concentration used for calculations and can be
#' configured to store raw absorptance values. However, because these are not
#' always stored, and calculation from raw absorptance causes small deviations
#' to the concentrations, it is useful to store uncorrected values. For H2O
#' these are concentrations, for CO2 span-corrected absorptance per unit
#' pressure. These values can then be used in recalculate gas mol fractions
#' at different oxygen levels. The function is called on data loading and not
#' intended to be exported.
#'
#' Valid factory and head calibration information is needed.
#'
#' @param tb a gas-exchange tibble with SysConst.Oxygen, Meas.H2Oa, Meas.CO2a,
#'   Meas.H2Or, Meas.CO2r, Status.Ts and Status.Tr
#'
#' @returns a gas-exchange tibble with Raw.H2Oa, Raw.H2Or, Raw.CO2aAbsP,
#'   Raw.CO2rAbsP
#' @noRd
calculate_raw <- function(tb) {
  #FIXME: I think we can move this to calcs
  if (!all(c("Meas.H2Oa", "Meas.H2Or") %in% names(tb))) {
    return(tb)
  }
  hc <- get_cal(tb$SysConst.UserCal, tb$SysConst.FactCal)

  psih <- psiH2O(hc, tb$SysConst.Oxygen)
  tb$Raw.H2Oa <- (tb$Meas.H2Oa / psih)@"mmol*mol^-1"
  tb$Raw.H2Or <- (tb$Meas.H2Or / psih)@"mmol*mol^-1"

  opsi_a <- psiCO2(hc, tb$SysConst.Oxygen, tb$Meas.H2Oa)
  # note units cannot correctly convert 1/degC to 1/K, so convert manually
  tb$Raw.CO2aAbsP <- frac2abs(hc, "co2", "a",
                               tb$Meas.CO2a / (tb$Status.Ts@"K" * opsi_a)) *
    opsi_a

  opsi_b <- psiCO2(hc, tb$SysConst.Oxygen, tb$Meas.H2Or)
  tb$Raw.CO2rAbsP <- frac2abs(hc, "co2", "b",
                               tb$Meas.CO2r / (tb$Status.Tr@"K" * opsi_b)) *
    opsi_b

  tb

}

#diff <1e-6 in my test. Seems fine
#FIXME: not yet used!
# ab <-function(Icr, Icro, Iwr, Iwro, tirga) {
#
#   # from factory for irga_b:
#   z <- 3.05e-5
#   Xwcr <- -0.000933 #actually called Xch
#   h2obzero <- 1.1035 #from datafile
#   #the manual has mistakes in Eq. 9-33:
#   ab <- 1 - (Iwr/Iwro + Xwcr * (1 - Icr / Icro)) *
#     (h2obzero + z * drop_units(tirga))
#
#   #from factory cal irga_b:
#   s3 <- 35.1; s2 <- -0.000223; s3 <- 5.38e-5
#   ab <- ab * (1 + (as.numeric(tirga) - 35.1) * (5.38e-5 + ab * -0.000223))
#   ab
#
# }
