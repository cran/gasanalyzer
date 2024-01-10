
#' Reads GFS-3000 text files and creates a tibble with gas-exchange data
#'
#' The text files stored by the GFS-3000 contain measured and calculated values
#' that are read by this function and formatted in a large tibble for use with
#' R. Note that no recalculation of the gas-exchange parameters is
#' performed, although it is possible to do that using [recalculate()] after
#' importing the data.
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
    file_encoding = "ISO-8859-1"

  headunit <- read.delim(file = filename, sep = delim, nrows = 2,
                         fileEncoding = file_encoding, header = FALSE,
                         colClasses = "character", row.names = NULL)

  header <- as.character(headunit[1, ])
  header_units <- as.character(headunit[2, ])

  # checks
  if (!all(c("H2Oabs", "dH2OZP", "dH2OMP", "Code") %in% header) ||
      !(("ppm") %in% header_units)) {
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

  df <- within (df, {
    Id <- seq.int(nrow(df))
    ETR <-  g0("ETR", NA_real_, "\U00B5mol*m^-2*s^-1")
    Imp <- g0("Imp", NA_real_, "steps")
    Tcuv <- g0("Tcuv", NA_real_, "degC")
    Tleaf <- g0("Tleaf", NA_real_, "degC")
    PARtop <- g0("PARtop", NA_real_, "\U00B5mol*m^-2*s^-1")
    Comment <- get0("Comment", ifnotfound = NA_character_)

    Aux1 <- g0("Aux1", 0, "mV")
    Aux2 <- g0("Aux2", 0, "mV")

    Oxygen <- gfs_calc_o2(get0("H2Obuf", ifnotfound = get0("H2Oabs")),
                                           get0("wa"), get0("dH2OMP"),
                                    get0("dH2OZP"))
  })

  testo2 <- vec_fill_missing(comment_to_oxygen(df$Comment))

  if (any(abs(testo2 - df$Oxygen) > 1, na.rm = T))
    warning("Oxygen concentration in comment field does not alway match the",
            "concentration back-calculated from the water mol fractions.\n")
  #trust comment, then back-calc
  df$Oxygen <- set_units(blend(testo2,df$Oxygen), "%")

  if (!unified_names) {
    names(df) <- make.names(names(df))
    as_tibble(df)
  } else {
    colnames(df) <- rename_header(names(df), "GFS3000")
    #deal with matching columns:
    splitcode <- stri_split_regex(df$GFS3000.Code, "(?<=[A-Za-z])(?=[0-9])|_",
                                  simplify = TRUE)
    # To make the data tidy, we need to restructure ZP(match) data
    # to match measurements to which it belongs
    zptrue <- grepl("ZP", splitcode[ , 1], fixed = TRUE)
    na_vec <- rep(NA, length(zptrue))
    o2fac <- gfs_o2_factor(df$Raw.H2Or, df$SysConst.Oxygen)

    #g0 instead of get0 is used in within because math fails if NULL w/o units
    #we also keep using ppm as units might (wrongly!) simply mmol/mol

    df <- within(df, {
      #GFS3000.Code <- splitcode[ , 1]
      SysObs.Averaging <- set_units(as.numeric(splitcode[ , 2]), "s")
      # why do none of the instruments save a timezone?
      SysObs.Time <- paste(get0("SysObs.Date"), get0("SysObs.HHMMSS")) |>
        as.POSIXct(tz = tz)
      SysObs.Instrument <- "GFS3000"
      SysObs.Filename <- bname
      GasEx.Ca <- get0("Meas.CO2s")
      GasEx.Time <- SysObs.Time - as.numeric(SysObs.Averaging) / 2

      Meas.H2Or <- o2fac * g0("Raw.H2Or", NA, "ppm")
      Meas.CO2a <- g0("Meas.CO2r", NA, "ppm") + g0("GFS3000.dCO2MP", NA, "ppm")
      Meas.H2Oa <- Meas.H2Or + g0("GFS3000.dH2OMP", NA, "ppm") * o2fac
      Raw.H2Oa <- Meas.H2Oa / o2fac
      # I think all GFS chambers have only one TC
      LTConst.fT1 <- 1
      LTConst.fT2 <- 0
      LTConst.fTEB <- 0
      # 0 sometimes means NA
      FLR.Fv_Fm <- get0("Fv_Fm", ifnotfound = NA_real_)
      FLR.Fvp_Fmp <- get0("Fvp_Fmp", ifnotfound = NA_real_)
      FLR.Fv_Fm[FLR.Fv_Fm == 0] <- NA_real_
      FLR.Fvp_Fmp[FLR.Fvp_Fmp == 0] <- NA_real_
      LeafQ.Qin <- g0("LeafQ.Qin", get0("Meas.QambIn", ifnotfound = NA),
                      "\U00B5mol*m^-2*s^-1")

      #Walz version of "Dynamic"
      if (exists("GFS3000.H2Obuf")) {
        Const.UseDynamic <- TRUE
        Dynamic.Hr <- o2fac * get0("GFS3000.H2Obuf")
        Meas.H2Os <- g0("Meas.H2Os", NA, "ppm")
        Meas.Flow <- g0("Meas.Flow", NA, "\U00B5mol*s^-1")
        Const.S <- g0("Const.S", NA, "cm^2")
        #might want to this from predefined eq instead:
        Dynamic.Edyn <- Meas.Flow * (Meas.H2Os - Dynamic.Hr) /
           (Const.S * (unity - Meas.H2Os))
        #CO2buf is rare, maybe only in some GFS versions?
        if (exists("Dynamic.Crd")) {
          Dynamic.Adyn <- (g0("Dynamic.Crd", NA, "ppm") -
                             g0("Meas.CO2s", NA, "ppm") *
             (unity - Dynamic.Hr) / (unity - Meas.H2Os)) * Meas.Flow /
            Const.S
        }
        rm("GFS3000.H2Obuf")
      } else
        Const.UseDynamic <- FALSE

      MchEvent.CO2match <- -g0("GFS3000.dCO2ZP", NA, "ppm")
      MchEvent.H2Omatch <- -g0("GFS3000.dH2OZP", NA, "ppm") * o2fac
      Meas2.dCO2 <- g0("GFS3000.dCO2MP", NA, "ppm") + MchEvent.CO2match
      Meas2.dH2O <- g0("GFS3000.dH2OMP", NA, "ppm") * o2fac + MchEvent.H2Omatch

      rm("GFS3000.dCO2ZP", "GFS3000.dH2OZP", "GFS3000.Code", "GFS3000.dCO2MP",
         "GFS3000.dH2OMP")
    })

    # Unfortunately, we have to fist initialize new cols
    df[c("MchEvent.Time", "MchEvent.HHMMSS", "MchEvent.CO2at",
         "MchEvent.H2Oat", "MchEvent.CO2adj", "MchEvent.H2Oadj",
         "MchEvent.Averaging", "MchStatus.Status")] <- list(NA)
    #and assign them only to zptrue rows
    df[zptrue, ] <- within(df[zptrue, ], {
      MchEvent.Time <- get0("SysObs.Time")
      MchEvent.HHMMSS <- get0("SysObs.HHMMSS")
      MchEvent.CO2at <- get0("Meas.CO2r")
      MchEvent.H2Oat <- get0("Meas.H2Or")
      MchEvent.CO2adj <- g0("Meas.CO2r", NA, "ppm") - g0("Meas.CO2s", NA, "ppm")
      MchEvent.H2Oadj <- g0("Meas.H2Or", NA, "ppm") -
        g0("Meas.H2Os", NA, "ppm")
      MchEvent.Averaging <- get0("SysObs.Averaging")
      MchStatus.Status <- get0("GFS3000.Status")
    })

    tmp <- grepl("MchEvent", colnames(df), fixed = T)
    df[tmp] <- filldown(df[tmp])
    df <- df[!zptrue, ] |> as_tibble()
    #reassign Obs and elapsed
    df[["SysObs.Obs"]] <- seq.int(nrow(df))
    df[["SysObs.Elapsed"]] <- as.numeric(df[["SysObs.Time"]] -
                                           df[["SysObs.Time"]][1])
    df <- fixup_import(df)

    units_options("simplify" = old_opt)
    df[sort_names(names(df))]
  } #end unified names
}
