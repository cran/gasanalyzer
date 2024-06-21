
# an environment for internal data
.gasanalyzerEnv <- new.env(parent = emptyenv())

#' Loads all internal data into a special environment
#'
#' @importFrom stats setNames
#' @noRd
load_internal_data <- function() {

  vars <- asNamespace("gasanalyzer")$vars
  #TODO: could consider moving the eqs assignment also to raw-data...
  eqs <- calcs()
  for (x in names(eqs)) {
    vars$fn[vars$fullname == x] <- list(eqs[x])
  }
  #no need to keep it around:
  rm("vars", envir = asNamespace("gasanalyzer"))

  # all these are fancy looking, but terrible to actually work with
  sympairs <- c(
    "\u0394" , "d",     # delta symbol
    "\u03b1" , "alpha", # alpha symbol
    "\'"     , "p",     # prime
    "\""     , "",      # unquote
    " "      , "",      # cumbersome
    "?"      , "",      # not needed
    "/"      , "_",     # silly
    "("      , "_",     #
    ")"      , "_",     #
    "@"      , "_at_",  #
    ":"      , "_",     # why do you torment us so?
    "-"      , "_",     #
    "%"      , "pct"    #
  )

  # A table to translate XL functions to R functions. For safety,
  # ONLY functions listed here will be parsed. So If new
  # XL function are used in a file, the need to be added here
  # or parsing formulas will be disabled!
  xl_to_r <- c("IF" = "ifelse", "LOG10" = "log10", "LN" = "log",
               "(?<![<>=])=" = "==", "<>" = "!=", "MAX" = "pmax",
               "MIN" = "pmin", "LEFT" = "left", "RIGHT" = "right",
               "SQRT" = "sqrt", "SIGN" = "sign", "EXP" = "exp",
               "ABS" = "abs", "COS" = "cos", "SIN" = "sin", "TAN" = "tan",
               "POWER" = "`^`", "PI" = "pi_xl", "CEILING" = "ceiling_xl",
               "DEGREES" = "degrees", "RADIANS" = "radians", "LOG" = "log")

  # useful parts from the calibration data
  needed_cals <- list(hc = c("co2aspan1", "co2aspan2", "co2bspan1",
                             "co2bspan2", "h2oaspan1", "h2oaspan2",
                             "h2obspan1", "h2obspan2"),
                      fc = c("bb.ch", "bb.cx", "bb.hx0", "bb.hx1",
                             "bb.pcorr_c1", "bb.pcorr_c2", "bb.pcorr_c3",
                             "bb.pcorr_h1", "bb.pcorr_h2", "bb.pcorr_h3",
                             "irga_b.Xch", "irga_a.co2.a1", "irga_a.co2.a2",
                             "irga_a.co2.a3", "irga_a.co2.a4", "irga_a.co2.a5",
                             "irga_a.co2.z", "irga_a.co2.s1", "irga_a.co2.s2",
                             "irga_a.co2.s3", "irga_b.co2.a1", "irga_b.co2.a2",
                             "irga_b.co2.a3", "irga_b.co2.a4", "irga_b.co2.a5",
                             "irga_b.co2.z", "irga_b.co2.s1", "irga_b.co2.s2",
                             "irga_b.co2.s3", "irga_a.h2o.a1", "irga_a.h2o.a2",
                             "irga_a.h2o.a3", "irga_a.h2o.z", "irga_a.h2o.s1",
                             "irga_a.h2o.s2", "irga_a.h2o.s3", "irga_b.h2o.a1",
                             "irga_b.h2o.a2", "irga_b.h2o.a3", "irga_b.h2o.z",
                             "irga_b.h2o.s1", "irga_b.h2o.s2", "irga_b.h2o.s3",
                             "cal.time", "caldate", "serialnumber"))

  # Assign all to a specific environment
  assign("xl_to_r", xl_to_r, envir = .gasanalyzerEnv)
  assign("sympairs", sympairs, envir = .gasanalyzerEnv)
  assign("needed_cals", needed_cals, envir = .gasanalyzerEnv)
  assign("vars", vars, envir = .gasanalyzerEnv)
}
