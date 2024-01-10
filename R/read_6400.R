#' Reads 6400XT text files and creates a tibble with gas-exchange data.
#'
#' The text files stored by the 6400 contain measured and calculated values that
#' are read by this function and formatted in a large tibble for use with R.
#' Constants and metadata are also added as columns. Note that no recalculation of
#' the gas-exchange parameters is performed, although it is possible to do that using
#' [recalculate()] after importing the data.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match.
#'
#' @param filename an text file containing 6400XT gas-exchange data.
#' @param tz a character string specifying the timezone for the loaded file. If
#'   omitted, the current time zone is used. Invalid values are typically
#'   treated as UTC, on some platforms with a warning.
#'
#' @returns A tibble with gas-exchange data in columns.
#'
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
#'   stri_split_fixed stri_detect_regex
#' @importFrom tibble tibble as_tibble
#' @importFrom units set_units units_options
#' @importFrom xml2 read_html as_list
#' @importFrom tools file_path_sans_ext
#' @importFrom utils compareVersion
#'
#' @seealso recalculate
#' @export
#'
#' @examples
#' example <- system.file("extdata//6400-testfile", package = "gasanalyzer")
#'
#' # read data
#' li6400data <- read_6400_txt(example)
#'
#' #View
#' li6400data
#'
read_6400_txt <- function(filename, tz = Sys.timezone()) {

  # there are various faster options, but readLines isn't that slow for
  # normal file sizes. Note that some files seem null-truncated...
  rawfile <- readLines(filename, skipNul = TRUE)

  if (!grepl("^\"OPEN ", rawfile[1]))
    stop(filename, " does not appear to be a 6400XT OPEN file.")

  # find sections. search only first 1000 lines for header, and
  # assume header < 1000 lines
  rfl <- length(rawfile)
  data_start <- which(rawfile[1:min(rfl, 1000)] == "$STARTOFDATA$")
  if (length(data_start) == 0L) {
    warning("No STARTOFDATA section in first 1000 lines, ignoring ",
            filename, ". Report a bug if this file is a valid OPEN file.\n")
    return(tibble())
  }

  tryCatch({
    # failing to read metadata might still allow reading data...
    # TODO: using cal data
    metadata <- read_html(paste0(rawfile[3:data_start - 2],
                                 collapse = "")) |> as_list()

    open_vers <- metadata[["html"]][["body"]][["p"]][["open"]][["version"]][[1]]
    open_vers <- gsub("\"", "", open_vers, fixed = TRUE)

    if (length(open_vers) == 0 ||
        is.na(numeric_version(open_vers, strict = FALSE)) ||
        compareVersion(open_vers, "6.2.0") < 0)
      warning(filename,
              " was created by an untested firmware version (v",
              open_vers, ").\n")
  },
  error = function(cond) {
    message("Failed to read metadata from ", filename, ":")
    message(cond)
    metadata <- NA
  })

  # lack of tz info is cumbersome
  filedate <- rawfile[2] |>
    gsub("\"","",x = _) |>
    gsub("Thr","Thu", x = _) |>
    as.POSIXct(tz = tz, format="%a %b %d %Y %H:%M:%S")

  rawdata <- rawfile[(data_start + 1):rfl]
  # data start with a number
  datarows <- stri_detect_regex(rawdata,"^[[:digit:]]+\t")
  #add header:
  datarows[1] = TRUE

  if (!any(datarows)) {
    warning("No data rows in ", filename, ".\n")
    return(tibble())
  }

  datamat <- stri_split_fixed(rawdata[datarows],
                              pattern = "\t", simplify = TRUE)

  # remove cols that are completely empty (mainly to work around bugs
  # happening during file creation...)
  datamat <- datamat[ , colSums(datamat != "") != 0]
  # remaining "" is really NA
  datamat[datamat == ""] <- NA

  header <- datamat[1, ] |>
    gsub("\"","", x = _) |>
    rename_header("Li6400") |>
    make.unique()
  # and apply:
  colnames(datamat) <- header

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  out <- datamat[-1,] |>
    as_tibble()

  out <- within(out, {
    SysObs.Instrument <- "6400"
    SysObs.Filename <-  file_path_sans_ext(basename(filename))
    SysObs.HHMMSS <- gsub("\"","", get0("SysObs.HHMMSS", ifnotfound = NA))
    FlrLS.fblue <- as.numeric(get0("Li6400.pctBlue", ifnotfound = NA)) / 100
    if (exists("Li6400.pctBlue")) rm("Li6400.pctBlue")
    SysObs.Time <- paste(as.Date(filedate), SysObs.HHMMSS) |>
      as.POSIXct(tz = tz)

    SysObs.Date <- as.character(as.Date(SysObs.Time))
    MchEvent.Time <- SysObs.Time - get0("Li6400.mchElpsd", ifnotfound = NA)
    #6400 doesn't store fan speeds:
    Meas.FanSpeed <- NA
    LTConst.fT1 <- 1.0 - as.numeric(get0("LTConst.fTEB", ifnotfound = "0"))
    LTConst.fT2 <- 0.0
    Meas.CO2s <- get0("Meas.CO2s", ifnotfound = NA)
    GasEx.Ca <- get0("GasEx.Ca", ifnotfound = Meas.CO2s)
    GasEx.E <- g0("GasEx.Emm", NA, "mol*m^-2*s^-1")
    LeafQ.Qin <- g0("LeafQ.Qin", get0("Meas.QambIn", ifnotfound = NA),
                    "\U00B5mol*m^-2*s^-1")
  })


  out <- fixup_import(out)

  units_options("simplify" = old_opt)
  out[sort_names(names(out))]

}
