#' Reads CIRAS-4 csv files and creates a tibble with gas-exchange data
#'
#' The csv files stored by the CIRAS-4 contain measured and calculated values
#' that are read by this function and formatted in a large tibble for use with
#' R. Note that no recalculation of derived variables are performed, although it
#' is possible to do so using [recalculate()] after importing the data.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match.
#'
#' @param filename a csv file containing gas-exchange data.
#'
#' @returns a tibble with gas-exchange data in columns.
#'
#' @importFrom stringi stri_list2matrix stri_trim stri_split_fixed
#' @importFrom tibble tibble as_tibble
#' @importFrom units set_units units_options
#' @importFrom tools file_path_sans_ext
#' @importFrom stats setNames
#' @seealso [recalculate()]
#'
#' @export
#' @examples
#' example <- system.file("extdata", "ciras4.csv", package = "gasanalyzer")
#'
#' # Read using unified column names:
#' cir4 <- read_ciras4(example)
#'
#' # Recalculate data using default gas exchange equations:
#' cir4_recalc <- recalculate(cir4, create_equations(c("default", "ciras4")))
#'
#' # View differences:
#' all.equal(cir4, cir4_recalc[names(cir4)], tol = 0.001)
read_ciras4 <- function(filename) {

  # for multiple files, we need to make sure every filename is unique
  # for now this is just a place-holder
  bname <- file_path_sans_ext(basename(filename))

  tryCatch( {
    # files seem asci, but actually pass utf8, i will not add an explicit
    # encoding for now
    ciras <- read.csv(filename)
  }, error = function(cond) {
    message("Failed to read ", filename, ":")
    message(cond, "\n", appendLF = F)
  })

  if (!exists("ciras")) {
    #warning("Failed to read ", filename, ".")
    return(tibble())
  }

  # checks
  if (length(nrow(ciras)) < 1 || nrow(ciras) < 1) {
    warning(filename, "contains no data.\n")
    return(tibble())
  }
  if (!all(c("TIMESTAMP", "DateTime", "H2Oa", "H2Or", "CO2a", "CO2r", "Aleaf",
             "Flow") %in% names(ciras))) {
    warning("Invalid file header, or no PLC4 cuvette attached, ",
            "ignoring ", filename, ".\n")
    return(tibble())
  }

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  #apply ciras units:
  nms <- names(ciras)
  mtch <- match(nms, .gasanalyzerEnv$vars$CIRAS4)
  mtch <- mtch[!is.na(mtch)]
  cunits <- setNames(object = .gasanalyzerEnv$vars$CIRAS4units[mtch],
                     .gasanalyzerEnv$vars$CIRAS4[mtch])
  ciras <- units_convert(ciras, cunits)

  # adjust columns that cannot be handled by the translation table:
  dt <- ciras[["DateTime"]] |>
    stri_trim() |>
    stri_split_fixed(pattern=" ") |>
    stri_list2matrix(byrow = T, n_min = 2)

  ciras["SysObs.Date"] <- dt[ ,1]
  ciras["SysObs.HHMMSS"] <- dt[ ,2]
  ciras["SysObs.Instrument"] <- "CIRAS4"
  ciras["SysObs.Filename"] <- bname
  ciras["Meas.FanSpeed"] <- NA

  #TEMPERATURE
  ciras[c("LTConst.fT1", "LTConst.fT2", "LTConst.fTEB")] <- list(1, 0, 0)
  ebrows <- stri_trim(ciras[["Tsensor"]]) == "EB"
  ciras[ebrows, "LTConst.fTEB"] <- 1
  ciras[ebrows, "LTConst.fT1"] <- 0

  #init light:
  ciras[c("QConst.fQambIn", "QConst.fQambOut", "QConst.fQflr",
          "QConst.fQheadLS")] <- 0

  ciras["Lcontrol"] <- g0("Lcontrol", NA_character_, envir = list2env(ciras))
  # we depend on this being in % later, so force it:
  ciras["FarRed"] <- g0("FarRed", 0, "%", TRUE, envir = list2env(ciras))

  ciras[c("HeadLS.Q", "Meas.QambIn", "FlrLS.Q")] <- list(
    set_units(0, "\U00B5mol*m^-2*s^-1"))

  #gbw in the package is defined as a one-sided conductance:
  ciras["Const.CustomBLC"] <- 1 / (2 * g0("rb", NA_real_, "m^2*s*mol^-1",
                                    envir = list2env(ciras)))
  # we have no idea how it depends on fan on area, so fixing gbw:
  ciras["GasEx.gbw"] <- ciras["Const.CustomBLC"]
  #Q incl. farred:
  PFD <- g0("PARi", NA_real_, "\U00B5mol*m^-2*s^-1",
                     envir = list2env(ciras)) /
    (1 - as.numeric(ciras[["FarRed"]]) / 100)


  # manual implies that the mass flow meter has been calibrated to this vol
  # at 0 degC and STP, therefore, since the mass flow scales linear with mol
  # flow, the molar volume at 0 and STP should be used:
  ciras["Meas.Flow"] <- set_units(g0("Flow", NA_real_, "cc*min^-1",
                                     envir = list2env(ciras)), "L*s^-1") /
    set_units(22.414, "\U00B5L*\U00B5mol^-1")



  #LIGHT
  # ciras defaults to 0.84 in manual and do no spectrum weighing. We therefore
  # use red for all calcs for now
  ciras["LQConst.RedAbsLED"] <- 0.84
  ciras["LQConst.RedAbsFlr"] <- 0.84
  # but i think this value is unlikely for ambient, nevertheless i stuck to it:
  ciras["LQConst.AbsAmbient"] <- 0.84

  # for LEDS the manual uses .14 including abs, and do not spectrum weighing.
  # I use red to reproduce this:
  ciras["LQConst.RedConvLED"] <- 0.14 / ciras["LQConst.RedAbsLED"]
  # the manual doesn't mention clearly what factors are used for the fluorometer
  # perhaps the same as for the LEDs..this would be about .16, which is
  # likely when the following is used:
  ciras["LQConst.RedConvFlr"] <- 0.14 / ciras["LQConst.RedAbsFlr"]
  # the manual claims 0.24 for ambient, incl abs, however this would require
  # abs=0.9. this is unlikely for ambient. 6800 is using .1911 which is quite
  # low the official value should be 0.217 or so....
  # perhaps there are some unexplained corrections, so i stuck to manual:
  ciras["LQConst.ConvAmbient"] <- 0.24 / ciras["LQConst.AbsAmbient"]

  color_cols <- c("Red", "Green", "Blue", "White", "FarRed")
  if (all("F" %in% names(ciras))) {
    ciras["FlrLS.Q"] <- PFD
    color_cols_new <- paste0("FlrLS.f", tolower(color_cols))
    ciras["QConst.fQflr"] <- 1
  } else { #FIXME: how about colors for ambient light?
    LEDrows <- stri_trim(ciras[["Lcontrol"]]) == "LED"
    ciras[LEDrows, "QConst.fQheadLS"] <- 1
    ciras[!LEDrows, "QConst.fQambIn"] <- 1
    ciras[LEDrows, "HeadLS.Q"] <- PFD[LEDrows]

    ciras[!LEDrows, "Meas.QambIn"] <- PFD[!LEDrows] *
      (1 - as.numeric(ciras[!LEDrows, "FarRed"]) / 100)
    color_cols_new <- paste0("HeadLS.f", tolower(color_cols))
  }

  # the NQP thing is probably a bug:
  ciras <- rename_cols(ciras,
                       c("H2Od", "CO2d", "NQP", color_cols),
                       c("Meas2.dH2O", "Meas2.dCO2", "NPQ", color_cols_new))
  ciras[c("DateTime", "Flow", "rb", "PARi")] <- list(NULL)

  colnames(ciras) <- rename_header(colnames(ciras), "CIRAS4")

  # final adjustments:
  ciras[["Meas.QambOut"]] <- g0("Meas.QambOut", 0, "\U00B5mol*m^-2*s^-1",
                                envir = list2env(ciras))
  ciras["LeafQ.Qin"] <-
    ciras[["QConst.fQambIn"]] * ciras[["Meas.QambIn"]] +
    ciras[["QConst.fQambOut"]] * ciras[["Meas.QambOut"]] +
    ciras[["QConst.fQflr"]] * ciras[["FlrLS.Q"]] *
    (set_units(1, "1") - g0("FlrLS.ffarred", 0, "1", envir = list2env(ciras))) +
    ciras[["QConst.fQheadLS"]] * ciras[["HeadLS.Q"]] *
    (set_units(1, "1") - g0("HeadLS.ffarred", 0, "1", envir = list2env(ciras)))

  #only light source weighted, not spectrum weighted:
  ciras["LeafQ.alpha"] <-
    (ciras[["QConst.fQambIn"]] * ciras[["LQConst.AbsAmbient"]] +
       ciras[["QConst.fQambOut"]] * ciras[["LQConst.AbsAmbient"]] +
       ciras[["QConst.fQflr"]] * ciras[["LQConst.RedAbsFlr"]] +
       ciras[["QConst.fQheadLS"]] * ciras[["LQConst.RedAbsLED"]]) /
    (ciras[["QConst.fQambIn"]] + ciras[["QConst.fQambOut"]] +
       ciras[["QConst.fQflr"]] + ciras[["QConst.fQheadLS"]])

  ciras["LeafQ.Conv"] <-
    (ciras[["QConst.fQambIn"]] * ciras[["LQConst.ConvAmbient"]] +
       ciras[["QConst.fQambOut"]] * ciras[["LQConst.ConvAmbient"]] +
       ciras[["QConst.fQflr"]] * ciras[["LQConst.RedConvFlr"]] +
       ciras[["QConst.fQheadLS"]] * ciras[["LQConst.RedConvLED"]]) /
    (ciras[["QConst.fQambIn"]] + ciras[["QConst.fQambOut"]] +
       ciras[["QConst.fQflr"]] + ciras[["QConst.fQheadLS"]])

  ciras["GasEx.Ca"] <- g0("Meas.CO2s", NA_real_, "\U00B5mol*mol^-1",
                          envir = list2env(ciras))
  ciras["GasEx.E"] <- g0("GasEx.Emm", NA, "mol*m^-2*s^-1",
                         envir = list2env(ciras))

  ciras <- ciras |> fixup_import() |> as_tibble()


  units_options("simplify" = old_opt)
  return(ciras[sort_names(names(ciras))])

}
