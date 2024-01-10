#' Write a gas-exchange tibble to a text file.
#'
#' The column names and column units are saved as a two-row header. The files
#' use UTF-16LE and CRLF line headings for compatibility. By default, tabs are
#' used as delimiter. If you intend to open the file in a spreadsheet program,
#' it may be helpful to use `csv` as file extension.
#'
#' Note that for data-exchange between R sessions, [saveRDS()] and [readRDS()]
#' are faster, and support saving and loading list columns (such as calibration
#' information and equations). This method is primarily meant for exchanging
#' data with other software packages.
#'
#' @param df a tibble with gas-exchange data
#' @param filename path to the output file
#' @param delim delimiter to use for the file
#'
#' @importFrom utils write.table
#' @returns No return value. If there a problems writing the file, a warning or
#'   error will be shown.
#'
#' @export
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' example <- system.file("extdata//d13C.tsv", package = "gasanalyzer")
#'
#' # read data and recalculate using default gas-exchange equations:
#' df <- read_gasexchange(example) |>
#'   recalculate(create_equations("default"))
#'
#' # write recaculated data
#' write_gasexchange(df, "d13C_recalculated.tsv")
#' \dontshow{setwd(.old_wd)}
write_gasexchange <- function(df, filename, delim = "\t") {

  #When spreadsheet support for JSON matures, we can add a json export
  #which will probably allow lists to be saved.
  if (!inherits(df, "data.frame") || nrow(df) < 1 ||
      !all(c("SysObs.Time", "Meas.CO2r", "Meas.CO2s", "Meas.H2Or",
             "Meas.H2Os", "Meas.Flow") %in% names(df)))
  {stop("df needs to be a data frame with 1 or more rows, and its columns ",
        "need to include: SysObs.Time, Meas.CO2r, Meas.CO2s, Meas.H2Or, ",
        "Meas.H2Os, and Meas.Flow")}

  unit_header <- lapply(df, get_units)

  # some classes need to be written differently
  df[unit_header == "posix"] <- sapply(df[unit_header == "posix"], as.numeric)

  # although lists can be converted to characters and written, I feel
  # it does not make sense at the moment because converting the equations
  # to char list is going to be a parsing disaster. RDS is much more useful
  # as a tool to save and load R objects.
  if ("list" %in% unit_header) {
    warning("List columns were omitted. Use saveRDS() and readRDS()",
            "to exchange data in R.\n")
    for (x in names(df[unit_header == "list"])) {
      df[[x]] <- NA
    }
  }

  #NB: base R gave me lots of issues writing utf16 with just one BOM
  #at the start, i'm abusing stri_write_lines instead. The downshot
  #is that it appends an extra tab at the end of every line, essentially
  #appending an empty column.
  fenc <- "UTF16-LE"
  fd <- file(filename, "wb", encoding = fenc)
  tryCatch({
    #manually write a BOM in an attempt to support broken software...
    writeBin(as.raw(c(0xff, 0xfe)), fd)
    stringi::stri_write_lines(names(df), encoding = fenc,
                              sep = delim, con = fd)
    writeBin(as.raw(c(0x0d, 0,0x0a, 0)), fd) #EOL
    stringi::stri_write_lines(unit_header, encoding = fenc,
                              sep = delim, con = fd)
    writeBin(as.raw(c(0x0d, 0,0x0a, 0)), fd)
    rubish <- by(df, seq_len(nrow(df)), \(row) {
      stringi::stri_write_lines(row, encoding = fenc,
                                sep = delim, con = fd)
      writeBin(as.raw(c(0x0d, 0,0x0a, 0)), fd)
    })
  }, finally = {close(fd)})

}
#' Read gas-exchange data from a text file.
#'
#' Data stored by [write_gasexchange()] can be read by this method. The first
#' row is the header, the second row specify the units. File encoding must be
#' UTF-16LE (use the export as unicode txt option in Microsoft Excel).
#'
#' @param filename path to the input file
#' @param delim delimiter to use for the file
#'
#' @importFrom utils read.delim
#' @importFrom tibble as_tibble
#' @returns a tibble with gas-exchange data
#'
#' @export
#' @examples
#'
#' example <- system.file("extdata//d13C.tsv", package = "gasanalyzer")
#'
#' # read data
#' read_gasexchange(example)
read_gasexchange <- function(filename, delim = "\t") {

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)
  # despite BOM, seems Mac needs LE here:
  fenc = "UTF-16LE"

  headunit <- read.delim(file = filename, sep = delim, nrows = 2,
                         fileEncoding = fenc, header = FALSE,
                         colClasses = "character", row.names = NULL)

  if (nrow(headunit) != 2L || ncol(headunit) < 6) {
    warning("Invalid header in ", filename, " returning an empty tibble.\n")
    return(tibble())
  }

  df <- read.delim(file = filename, sep = delim, skip = 2,
                   fileEncoding = fenc, header = FALSE,
                   colClasses = "character", row.names = NULL,
                   na.strings = c("", "NA"))

  if (length(nrow(df)) < 1 || nrow(df) == 0 || ncol(df) != ncol(headunit)) {
    warning(filename, " contains no valid data. Returning an empty tibble.\n")
    return(tibble())
  }

  # avoid appending an empty column
  if (all(headunit[ncol(headunit)] == "") && all(is.na(df[ncol(df)]))) {
    df <- df[-ncol(df)]
    headunit <- headunit[-ncol(headunit)]
  }

  header <- as.character(headunit[1, ])
  header_units <- as.character(headunit[2, ])
  if (!all(c("SysObs.Time", "Meas.CO2r", "Meas.CO2s",
             "Meas.H2Or", "Meas.H2Os", "Meas.Flow") %in% header)) {
    warning("Invalid header in ", filename, " returning an empty tibble.\n")
    return(tibble())
  }
  names(df) <- header
  names(header_units) <- header

  df <- df |> as_tibble(.name_repair = "unique") |>
    units_convert(header_units) |> fixup_import()

  units_options("simplify" = old_opt)
  df[sort_names(names(df))]

}

#' Export a subset of the data into the ESS-DIVE reporting format for leaf-level
#' gas exchange data.
#'
#' Ely et al. (2021) proposed a standardized nomenclature for reporting
#' gas-exchange data and metadata within the framework of the Environmental
#' System Science Data Infrastructure for a Virtual Ecosystem (ESS-DIVE)
#' repository. This method converts data frames or tibbles created with
#' `gasanalyzer` to this standardized format. Note that the scope of the
#' proposed standard is limited, and therefore only a subset of the data
#' is exported. Users should add relevant additional columns and provide
#' relevant metadata.
#'
#' If a filename is given as argument, the data is written into a comma
#' separated, UTF-8 encoded file without BOM and with CRLF line headings. In
#' addition, a data dictionary file ("_dd" is inserted at the end of the
#' filename, before the file extension). If no filename is provided, the
#' converted data is returned.
#'
#' @param df a tibble with gas-exchange data.
#' @param filename path to the output file. If none provided the method returns
#'   the data as a tibble.
#' @param extra_cols a character vector specifying additional columns (not
#'   specified in the standard specified by Ely et al. 2021) to include in the
#'   returned data or saved files.
#'
#' @returns Nothing if a filename is provided. Otherwise, a tibble with
#'   variables and headings specified by the ESS-DIVE gas-exchange standard
#'   is returned.
#'
#' @references Ely KS, Rogers A, Agarwal DA, et al (2021) *A reporting format
#'   for leaf-level gas exchange data and metadata.* Ecol Inform 61:101232.
#'   https://doi.org/10.1016/j.ecoinf.2021.101232
#'
#' @importFrom utils write.table write.csv
#' @importFrom stats setNames
#'
#' @export
#' @examples
#' example <- system.file("extdata//d13C.tsv", package = "gasanalyzer")
#'
#' # read data and recalculate using default gas-exchange equations:
#' df <- read_gasexchange(example) |>
#'   recalculate(create_equations("default"))
#'
#' # view df in ess_dive format:
#' export_ess_dive(df)
#'
#' \dontshow{.old_wd <- setwd(tempdir())}
#' # save the data and a data dictionary:
#' export_ess_dive(df, "ess_dive_test.csv")
#' # read and show the dictionary:
#' readLines("ess_dive_test_dd.csv")
#' \dontshow{setwd(.old_wd)}
export_ess_dive <- function(df, filename = "", extra_cols = NULL) {

  # Ely at al. mentions these need to be included:
  required_vars <- c("GasEx.A", "GasEx.Ci", "Meas.CO2s", "GasEx.gsw",
  "Meas.Pa", "LeafQ.Qin", "GasEx.RHcham", "GasEx.TleafCnd", extra_cols)

  if (!inherits(df, "data.frame") || nrow(df) < 1 ||
      !all(required_vars %in% names(df)))
  stop("df needs to be a data frame with 1 or more rows, and its columns ",
        "need to include: ", paste(required_vars, collapse = ", "))

  #get relevant translation table:
  vars <- .gasanalyzerEnv$vars[!is.na(.gasanalyzerEnv$vars[["ESSDIVE"]]), ]
  nms <- names(df)
  nms <- nms[nms %in% c(vars[["fullname"]], extra_cols)]

  idx <- vars$fullname %in% nms
  idx2 <- nms %in% vars[["fullname"]]
  newnames <- descriptions <- nms
  newnames[idx2] <- vars[["ESSDIVE"]][idx]
  descriptions[nms %in% .gasanalyzerEnv$vars$fullname] <-
    .gasanalyzerEnv$vars[.gasanalyzerEnv$vars$fullname %in% nms, "description"]

  #only enforce units specified by ESSDIVE, leave others as-is.
  #FIXME: c("HH:MM:SS", "YYYYMMDD") are not really units, and so
  #not known by units_convert. maybe i'll change that, for now
  #we'll apply these later.
  unitnames <- setNames(vars[["ESSDIVEunits"]][idx], newnames[idx2])

  out <- setNames(df[nms], newnames) |>
    units_convert(unitnames)

  #time and date for ESS-DIVE are always UTC, so use unixtime to set them:
  out[["time"]] <- format(df[["SysObs.Time"]], tz = "UTC", format = "%H:%M:%S")
  out[["date"]] <- format(df[["SysObs.Time"]], tz = "UTC", format = "%Y%m%d")
  #get units for all exported cols:
  unitnames <- sapply(out, get_units)
  #and set time and date the way ESSDIVE wants it:
  unitnames[c("time", "date")] <- c("HH:MM:SS", "YYYYMMDD")

  #missing values are N/A for text, but -9999 (for some reason) for numbers
  navals <- setNames(character(length(newnames)), newnames)
  for (x in newnames) {
    out[[x]][is.na(out[[x]])] <- navals[x] <-
      switch(class(out[[x]])[1], "numeric" = , "double"  = , "units"   = ,
             "integer" = -9999, "character" = "N/A", NA)
  }

  if (length(filename) == 1 && filename != "") {
    res <- try(write.csv(out, file = filename, quote = TRUE, eol = "\r\n",
                         row.names = FALSE, fileEncoding = "UTF8", na = "NA"))
    if (inherits(res,"try-error"))
      stop("Error writing data.")
    #also write a dict:
    hdr <- data.frame(V1 = c("Number_Header_Rows", "ColumnName_Row_Number",
                      "Column_Delimiter", "Datafile_ColumnName"),
                      V2 = c("1", "1", ",", "Datafile_Measure"),
                      V3 = c("", "", "", "Datafile_Units"),
                      V4 = c("", "", "", "Datafile_MissingValue"))

    unitnames[is.na(unitnames)] <- "N/A"
    data <- data.frame(V1 = newnames, V2 = descriptions,
                       V3 = unitnames, V4 = navals)
    # add dd to the filename, before the dot, or at the end:
    ddfilename <- gsub("^[^.]*?\\K([^.])(?=\\.|$)","\\1_dd",
                       filename, perl = TRUE)
    dd <- rbind(hdr, data)
    #we manually insert qoutes (write.table will otherwise quote all char vecs):
    for(x in names(dd)) {
      dd[[x]] <- gsub('"', '""', dd[[x]], fixed = TRUE)
      ndq <- grep('[\",]', dd[[x]])
      dd[ndq, x] <- paste0("\"", dd[ndq, x], "\"")
    }

    res <- try(write.table(dd, file = ddfilename, quote = FALSE,
                         eol = "\r\n", row.names = FALSE, fileEncoding = "UTF8",
                         na = "N/A", sep = ",", col.names = FALSE))

    if (!inherits(res,"try-error")) {
      message("Successfully written data to ", filename, " and dictionary to ",
              ddfilename, ".")
    } else
      stop("Error when writing dictionary ", filename, ".")
  }
  else
    out

}
