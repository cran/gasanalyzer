#' Reads CIRAS-3 xml files and creates a tibble with gas-exchange data
#'
#' The xml files stored by the CIRAS-3 contain a limited set of measured
#' data that are read by this function and formatted in a large tibble for use
#' with R. Note that no recalculation of the gas-exchange parameters is
#' performed. It is possible to calculate derived gas-exchange parameters
#' using [recalculate()], but required metadata (leaf area, etc) has to
#' be provided first.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match.
#'
#' @param filename an xml file containing gas-exchange data.
#' @param tz a character string specifying the timezone for the loaded file. If
#'   omitted, the current time zone is used. Invalid values are typically
#'   treated as UTC, on some platforms with a warning.
#'
#' @return a tibble with gas-exchange data in columns.
#'
#' @importFrom xml2 read_xml xml_children xml_attrs xml_name
#' @importFrom stringi stri_list2matrix
#' @importFrom tibble tibble as_tibble
#' @importFrom units set_units units_options
#' @importFrom tools file_path_sans_ext
#' @seealso recalculate
#'
#' @noRd
# TODO: I have no good ciras3 testfiles at the moment, export as soon as I have
# been able to test
read_ciras3_xml <- function(filename, tz = Sys.timezone()) {

  # always keep this off:
  units_options(simplify = NA)
  # for multiple files, we need to make sure every filename is unique
  # for now this is just a place-holder
  bname <- file_path_sans_ext(basename(filename))

  tryCatch( {
    ciras <- read_xml(filename)
    if (!isTRUE(xml_name(ciras) == "PpSystemsCiras3")) {
      warning("Not a Ciras3 datafile, ignoring ", filename, ".")
      return(tibble())
    }
    alldata <- xml_children(ciras) |>
      xml_attrs()
    data_length <- length(alldata)
    if (data_length < 1) {
      warning("No data in xml, ignoring ", filename, ".")
      return(tibble())
    }
  }, error = function(cond) {
    message("Failed to read ", filename, ":")
    message(cond)
    return(tibble())
  })
  #longest header may not suffice....:
  longest_header <- which.max(lengths(alldata))
  header <- names(alldata[[longest_header]]) |>
    rename_header("Ciras3")

  dt_mat <- stri_list2matrix(alldata, by_row = T)
  col_n <- ncol(dt_mat)

  if (col_n < 1) {
    warning("No data in xml, ignoring ", filename, ".")
    return(tibble())
  }

  if (!isTRUE(length(header) == col_n))
    warning("Header of ", filename, " may be incorrect.")

  colnames(dt_mat) <- header[1:col_n]
  df <- as_tibble(dt_mat, .name_repair="unique")

  mrows <- df[["Ciras3.RecType"]] == "M"
  if (length(mrows) != data_length) {
      warning("No RecType tag in xml? Ignoring ", filename, ".")
      return(tibble())
  }

  df <- df[mrows, ]
  row_n <- nrow(df)

  df <- within(df, {
    SysObs.obs <- seq(1, row_n)
    # default units are kPa but ciras3 reports mbar:
    Meas.Pa <- as.numeric(get0("Meas.Pa", ifnotfound = NA)) / 10
    Sysobs.Time <- as.POSIXct(as.numeric(get0("Ciras3.ExcelTime",
                                              ifnotfound = NA)) *
                                60 * 60 * 24, origin = "1899-12-30",
                              tz = "UTC") |> format() |>
      as.POSIXct(tz = tz)

    rm("Ciras3.RecType", "Ciras3.ExcelTime")
  }) |> fixup_import()

  return(df)

}
