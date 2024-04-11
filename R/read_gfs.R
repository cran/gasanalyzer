#' Reads GFS-3000 text files and creates a tibble with gas-exchange data
#'
#' The text files stored by the GFS-3000 contain measured and calculated values
#' that are read by this function and formatted in a large tibble for use with
#' R. Note that no recalculation of derived variables is performed, although it
#' is possible to do so using [recalculate()] after importing the data.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match.
#'
#' @param filename an xlsx file containing 6800 gas-exchange data.
#' @param tz a character string specifying the timezone for the loaded file. If
#'   omitted, the current time zone is used. Invalid values are typically
#'   treated as UTC, on some platforms with a warning.
#' @param unified_names = TRUE, use unified column names. This is necessary for
#'   further processing of the data using this package.
#' @param skip_to_data use skip=4 if the file has a double header.
#' @param delim = ";" Allows specified the delimiter used in the files. Re-saved
#'  data may use a comma as delimiter.
#'
#' @returns a tibble with gas-exchange data in columns and equations as
#'   attribute.
#'
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
#'   stri_split_fixed stri_detect_regex stri_enc_detect stri_split_regex
#' @importFrom tibble tibble as_tibble
#' @importFrom units set_units units_options
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.delim
#' @importFrom vctrs vec_fill_missing
#' @seealso [recalculate()]
#'
#' @export
#' @examples
#' example <- system.file("extdata//aci1.csv", package = "gasanalyzer")
#'
#' # Read using GFS-3000 names and formatting:
#' gfs3000_old <- read_gfs(example, unified_names = FALSE)
#' # Read using unified column names:
#' gfs3000 <- read_gfs(example)
#'
#' # Inspect the intercellular CO2:
#' gfs3000_old$ci
#' gfs3000$GasEx.Ci
#'
#' # Recalculate data using default gas exchange equations:
#' gfs3000 <- recalculate(gfs3000, create_equations(c("default", "gfs3000")))
#' gfs3000$GasEx.Ci
read_gfs <- function(filename, tz = Sys.timezone(), unified_names = TRUE,
                     skip_to_data = 2, delim = ";") {

  bname <- file_path_sans_ext(basename(filename))

  #guess file encoding:
  file_encoding <- stri_enc_detect(readChar(filename, 512))[[1]]$Encoding[1]

  if (is.na(file_encoding))
    file_encoding <- "ISO-8859-1"

  headunit <- read.delim(file = filename, sep = delim, nrows = 2,
                         fileEncoding = file_encoding, header = FALSE,
                         colClasses = "character", row.names = NULL)

  header <- as.character(headunit[1, ])
  header_units <- as.character(headunit[2, ])

  # checks
  if (!all(c("Date", "Time", "H2Oabs", "dH2OZP", "dH2OMP",
             "Code") %in% header) || !(("ppm") %in% header_units)) {
    warning("Invalid header, ignoring", filename, ".\n")
    return(tibble())
  }

  df <- read.delim(file = filename, sep = delim, skip = skip_to_data,
                   fileEncoding = file_encoding, header = FALSE,
                   colClasses = "character", row.names = NULL,
                   col.names = header, check.names = FALSE,
                   na.strings = c("", "NA"))
  # I think this reports errors:
  df <- df[df$Code != "BC", ]
  if (length(nrow(df)) < 1 || nrow(df) == 0) {
    warning(filename, "contains no data.\n")
    return(tibble())
  }

  names(header_units) <- header
  header_units[header %in% c("Date", "Time", "Object", "Code", "Status")] <- ""
  header_units[header_units == "-" | header_units == "mV"] <- ""

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  df <- units_convert(df, header_units)

  df["Id"] <- seq.int(nrow(df))
  df["ETR"] <-  g0("ETR", NA_real_, "\U00B5mol*m^-2*s^-1", envir = list2env(df))
  df["Imp"] <- g0("Imp", NA_real_, "steps", envir = list2env(df))
  df["Tcuv"] <- g0("Tcuv", NA_real_, "degC", envir = list2env(df))
  df["Tleaf"] <- g0("Tleaf", NA_real_, "degC", envir = list2env(df))
  df["PARtop"] <- g0("PARtop", NA_real_, "\U00B5mol*m^-2*s^-1",
                     envir = list2env(df))
  df["Comment"] <- g0("Comment", NA_character_, envir = list2env(df))

  df["Aux1"] <- g0("Aux1", 0, "mV", envir = list2env(df))
  df["Aux2"] <- g0("Aux2", 0, "mV", envir = list2env(df))

  df["Oxygen"] <- gfs_calc_o2(g0("H2Obuf", g0("H2Oabs", NA_real_,
                                              envir = list2env(df))),
                              g0("wa", NA_real_, envir = list2env(df)),
                              g0("dH2OMP", NA_real_, envir = list2env(df)),
                              g0("dH2OZP", NA_real_, envir = list2env(df)))

  testo2 <- vec_fill_missing(comment_to_oxygen(df[["Comment"]]))

  if (any(abs(testo2 - df[["Oxygen"]]) > 1, na.rm = T))
    warning("Oxygen concentration in comment field does not alway match the",
            "concentration back-calculated from the water mol fractions.\n")
  #trust comment, then back-calc
  df["Oxygen"] <- set_units(blend(testo2, df[["Oxygen"]]), "%")

  if (!unified_names) {
    names(df) <- make.names(names(df))
    as_tibble(df)
  } else {
    colnames(df) <- rename_header(names(df), "GFS3000")
    #deal with ZP / match columns:
    splitcode <- stri_split_regex(df[["GFS3000.Code"]],
                                  "(?<=[A-Za-z])(?=[0-9])|_",
                                  simplify = TRUE)
    # To make the data tidy, we need to restructure ZP(match) data
    # to match measurements to which it belongs
    zptrue <- grepl("ZP", splitcode[ , 1], fixed = TRUE)
    na_vec <- rep(NA, length(zptrue))
    o2fac <- gfs_o2_factor(df[["Raw.H2Or"]], df[["SysConst.Oxygen"]])

    df["SysObs.Averaging"] <- set_units(as.numeric(splitcode[ , 2]), "s")
      # why do none of the instruments save a timezone?
    df["SysObs.Time"] <- paste(df[["SysObs.Date"]], df[["SysObs.HHMMSS"]]) |>
        as.POSIXct(tz = tz)
    df["SysObs.Instrument"] <- "GFS3000"
    df["SysObs.Filename"] <- bname
    df["GasEx.Time"] <- df[["SysObs.Time"]] -
      as.numeric(df[["SysObs.Averaging"]]) / 2

    dfenv <- list2env(df)

    df["Meas.H2Os"] <- g0("Meas.H2Os", NA_real_, "mmol*mol^-1",
                          envir = dfenv)
    df["Raw.H2Or"] <- g0("Raw.H2Or", NA_real_, "mmol*mol^-1",
                         envir = dfenv)
    df["Meas.CO2r"] <- g0("Meas.CO2r", NA_real_, "\U00B5mol*mol^-1",
                          envir = dfenv)
    df["Meas.CO2s"] <- g0("Meas.CO2s", NA_real_, "\U00B5mol*mol^-1",
                         envir = dfenv)
    df["Meas.Flow"] <- g0("Meas.Flow", NA_real_, "\U00B5mol*s^-1",
                          envir = dfenv)
    df["Const.S"] <- g0("Const.S", NA_real_, "cm^2", envir = dfenv)
    df["Status.Status"] <- g0("Status.Status", NA_character_, envir = dfenv)
    df["FLR.Fv_Fm"] <- g0("FLR.Fv_Fm", NA_real_, envir = dfenv)
    df["FLR.Fvp_Fmp"] <- g0("FLR.Fvp_Fmp", NA_real_, envir = dfenv)
    df["GFS3000.dCO2MP"] <- g0("GFS3000.dCO2MP", NA_real_,
                               "\U00B5mol*mol^-1", envir = dfenv)
    df["GFS3000.dH2OMP"] <- g0("GFS3000.dH2OMP", NA_real_, envir = dfenv)
    df["MchEvent.CO2match"] <- -g0("GFS3000.dCO2ZP", NA_real_,
                                   "\U00B5mol*mol^-1", envir = dfenv)
    df["MchEvent.H2Omatch"] <- set_units(-g0("GFS3000.dH2OZP", NA_real_,
                                             envir = dfenv) * o2fac,
                                         "\U00B5mol*mol^-1")
    df["LeafQ.Qin"] <- g0("LeafQ.Qin", g0("Meas.QambIn", NA_real_,
                                          envir = dfenv),
                          "\U00B5mol*m^-2*s^-1", envir = dfenv)


    df["Meas.H2Or"] <- set_units(o2fac * df[["Raw.H2Or"]], "\U00B5mol*mol^-1")
    df["Meas.CO2a"] <- df[["Meas.CO2r"]] + df[["GFS3000.dCO2MP"]]
    df["Meas.H2Oa"] <- set_units(df[["Meas.H2Or"]] + df[["GFS3000.dH2OMP"]]
                                 * o2fac, "\U00B5mol*mol^-1")
    df["Raw.H2Oa"] <- set_units(df[["Meas.H2Oa"]] / o2fac, "\U00B5mol*mol^-1")
    df["Meas2.dCO2"] <- df[["GFS3000.dCO2MP"]] + df[["MchEvent.CO2match"]]
    df["Meas2.dH2O"] <- set_units(df[["GFS3000.dH2OMP"]] * o2fac +
                                    df[["MchEvent.H2Omatch"]], "mmol*mol^-1")

    df["GasEx.Ca"] <- df[["Meas.CO2s"]]

    # I think all GFS chambers have only one TC
    df[c("LTConst.fT1","LTConst.fT2", "LTConst.fTEB")] <- list(1, 0, 0)

   #Walz version of "Dynamic"
   if (all("GFS3000.H2Obuf" %in% names(df))) {
     df["Const.UseDynamic"] <- TRUE
     df["Dynamic.Hr"] <- o2fac * g0("GFS3000.H2Obuf", NA_real_,
                                    envir = list2env(df))

     df["Dynamic.Edyn"] <- df[["Meas.Flow"]] * (df[["Meas.H2Os"]] -
                                                      df[["Dynamic.Hr"]]) /
       (df[["Const.S"]] * (unity - df[["Meas.H2Os"]]))
      #CO2buf is rare, maybe only in some GFS versions?
      if (all("Dynamic.Crd" %in% names(df))) {
        df["Dynamic.Adyn"] <- (set_units(df[["Dynamic.Crd"]],
                                        "\U00B5mol*mol^-1") -
          df[["Meas.CO2s"]] * (unity - df[["Dynamic.Hr"]]) /
          (unity - df[["Meas.H2Os"]])) * df[["Meas.Flow"]] / df[["Const.S"]]
      }
      df["GFS3000.H2Obuf"] <- NULL
    } else
      df["Const.UseDynamic"] <- FALSE

    df[c("GFS3000.dCO2ZP", "GFS3000.dH2OZP", "GFS3000.Code", "GFS3000.dCO2MP",
         "GFS3000.dH2OMP")] <- NULL

    df[zptrue, "MchEvent.Time"] <- df[zptrue, "SysObs.Time"]
    df[zptrue, "MchEvent.HHMMSS"] <- df[zptrue, "SysObs.HHMMSS"]
    df[zptrue, "MchEvent.CO2at"] <- df[zptrue, "Meas.CO2r"]
    df[zptrue, "MchEvent.H2Oat"] <- df[zptrue, "Meas.H2Or"]
    df[zptrue, "MchEvent.CO2adj"] <- df[zptrue, "Meas.CO2r"] -
      df[zptrue, "Meas.CO2s"]
    df[zptrue, "MchEvent.H2Oadj"] <- df[zptrue, "Meas.H2Or"] -
      df[zptrue, "Meas.H2Os"]
    df[zptrue, "MchEvent.Averaging"] <- df[zptrue, "SysObs.Averaging"]
    df[zptrue, "MchStatus.Status"] <- df[zptrue, "Status.Status"]

    tmp <- grepl("MchEvent", colnames(df), fixed = T)
    df[tmp] <- filldown(df[tmp])
    df <- df[!zptrue, ] |> as_tibble()
    #reassign Obs and elapsed
    df[["SysObs.Obs"]] <- seq.int(nrow(df))
    df[["SysObs.Elapsed"]] <- as.numeric(df[["SysObs.Time"]] -
                                           df[["SysObs.Time"]][1])
    df <- fixup_import(df)
    # 0 sometimes means NA
    #df[!is.na(df[["FLR.Fv_Fm"]]) & df[["FLR.Fv_Fm"]] == 0, "FLR.Fv_Fm"] <-
    #  NA_real_
    #df[!is.na(df[["FLR.Fvp_Fmp"]]) & df[["FLR.Fvp_Fmp"]] == 0,"FLR.Fvp_Fmp"] <-
    #  NA_real_

    units_options("simplify" = old_opt)
    df[sort_names(names(df))]
  } #end unified names
}
