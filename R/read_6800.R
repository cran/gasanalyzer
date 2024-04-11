#' Read gas-exchange equations directly from 6800 xlsx files.
#'
#' It is recommended to use [create_equations()] for recalculating gas-exchange
#' data, but under some conditions, it may be useful to apply exactly the same
#' equations as used in the xlsx data file.
#'
#' Currently, this only works for xlsx files stored by the 6800. this function
#' extracts xlsx formulas from the file and stores them in a list for use by the
#' [recalculate()] function. Note there is no guarantee that the extracted
#' equations work on any other data files. Since newer versions of the 6800
#' firmware allows defining custom equations, it is not guaranteed that all
#' equations can be extracted successfully.
#'
#' @param filename an xlsx file containing 6800 gas-exchange data
#'
#' @returns A list of gas-equations.
#'
#' In principle, this can be made to work for the 6400 as well, but since
#' that instrument uses a variation of the older xls format, it is hard to
#' get working in practice.
#'
#' @seealso [create_equations()]
#' @export
#'
#' @examples
#' example <- system.file("extdata//lowo2.xlsx", package = "gasanalyzer")
#'
#' # get equations stored in the xlsx file
#' Eqs <- read_6800_equations(example)
#'
#' #Inpect how stomatal conductance is calculated:
#' Eqs$GasEx.gsw
read_6800_equations <- function(filename) {

  parse_6800_xlsx(filename, extract_formula = TRUE,
                  extract_units = FALSE)[["formula"]] |>
    extract_6800_equations()

}

#' Extract gas-exchange equations from parsed xlsx files
#'
#' Intended for internal use.
#'
#' @param df containing formulas, extracted from an xlsx using
#'   [parse_6800_xlsx()]
#'
#' @returns a list of quosures with equations.
#'
#' @importFrom utils getParseData
#' @importFrom stats setNames
#' @seealso [create_equations()]
#' @noRd
extract_6800_equations <- function(df) {

  # parsing input maybe unsafe. We check whether the functions are anticipated
  # and parsing was successful. Additional funcs can be added using xl_to_r
  tryCatch(
    { # use base to parse as we need the keep source:
      parsedForms <- parse(text = df$fn, keep.source = TRUE)
      pd <- getParseData(parsedForms)
      # all function calls have to be members of xl_to_r:
      fdif <- setdiff(pd$text[pd$token == "SYMBOL_FUNCTION_CALL"],
                      .gasanalyzerEnv$xl_to_r)
      if (!(length(fdif) == 0))
        stop("Please report unknown function calls: ", fdif, "\n")

      out <- setNames(as.list(parsedForms), nm = df$header)
      # add an entry to indicate these quos use no units:
      out[["gasanalyzer.UseEqUnits"]] <- FALSE
      return(out)
    },
    error = function(cond) {
      message("Error extracting XL formulas:")
      message(cond)
    }
  )
}

#' Parse 6800 xlsx files and return remarks, headers, units, formulas and data.
#'
#' Intended for internal use.
#'
#' Currently, this only works for xlsx files stored by the 6800. In principle,
#' it can be made to work for the 6400 as well, but since that instrument uses a
#' variation of the older xls format, it is hard to get working in practice.
#'
#' @param filename an xlsx file containing 6800 gas-exchange data
#'
#' @returns a list of dataframes with remarks, headers, units, formulas and data
#'
#' @importFrom tidyxl xlsx_cells
#' @importFrom utils getParseData modifyList
#' @importFrom stats setNames reshape
#' @importFrom vctrs vec_split vec_rbind vec_cbind
#' @importFrom utils compareVersion
#' @noRd
parse_6800_xlsx <- function(filename, extract_formula = FALSE,
                            extract_units = TRUE) {


  allCells <- xlsx_cells(filename, include_blank_cells = FALSE)
  class(allCells) <- "data.frame"

  if (any(allCells$sheet=="Measurements"))
    msheet <- "Measurements"
  else
    msheet <- allCells$sheet[1]

  # load head calibration from remarks
  remarksDF <- allCells[allCells$sheet == "Remarks" & allCells$col < 3L, ]
  remarksDF$X <- blend(remarksDF$character, as.character(remarksDF$numeric),
                       as.character(remarksDF$logical))
  remarksDF <- reshape(remarksDF[,c("row","col", "X")], direction="wide",
                       timevar = "col", idvar = "row", sep="", v.names="X")

  bluestem <- remarksDF$X2[remarksDF$X1 == "Console ver"]

  if (!grepl("Bluestem", bluestem, fixed = TRUE)) {
    warning(filename, " is not a valid 6800 xlsx file.\n")
    list(remarks = NA, headers = NA, units = NA,
         formula = NA, consts = NA, data = NA)
  } else {
    vers <- strsplit(bluestem,"v.")[[1]][2]
    if (is.na(numeric_version(vers, strict = FALSE)) ||
        compareVersion(vers, "2") < 0)
      warning(filename, " created by an unsupported firmware version (v",
              vers,").\n")
  }

  allCells <- allCells[allCells$sheet == msheet,
                       c("address", "row", "col", "error", "logical",
                         "character", "numeric", "formula", "data_type",
                         "content")]

  # this row has the header after which the measurement data is found:
  datastart <- allCells$row[!is.na(allCells$character) &
                              tolower(allCells$character) == "obs"] - 1L

  if (length(datastart) < 1L) {
    warning("No data rows in ", filename, ".")
    return(list(remarks = remarksDF, headers = NA, units = NA,
                formula = NA, consts = NA, data = NA))
  } else if (length(datastart) > 1L) {
    warning("Parse error in ", filename,
            ": Multiple rows with obs header unsupported! Ignoring data.\n")
    return(list(remarks = remarksDF, headers = NA, units = NA,
                formula = NA, consts = NA, data = NA))
  }

  datastart <- datastart[1L]

  # Measurement rows are identified by col 1 having a numeric value
  datarows <- allCells$row[!is.na(allCells$numeric) & allCells$col == 1L]
  headerrows <- rep(datastart, each = 3L) + 0:2L

  # create a df for our constants with names and addresses
  constsDA <- allCells[!(allCells$row %in% c(headerrows, datarows)), ]
  # consts blocks start on lines that have characters in col1
  constsDA$n <- cumsum(!is.na(constsDA$character) & constsDA$col == 1L)
  # split the df by grp, add col with the header name for the header row
  # and remove the original line with the header name
  constsDA <- lapply(vec_split(constsDA, constsDA$n)[[2]], function(x)
  { x$group <- replace(NA, x$row == x$row[1L],
                      x$character[1L]); x[-1L, ]}) |>
    do.call(vec_rbind, args = _)

  constsDF <- constsDA[!is.na(constsDA$group),
                       c("character", "group", "row", "col")]
  # create header names by merging the beheaded header
  # with the character cells to the right of them
  constsDF$header <- do.call(paste, c(constsDF[c("group", "character")],
                                      sep = ".")) |>
    rename_header("Li6800")

  # assumption: header is 1 row, data starts after that
  # and *is also only one row* (keep row integer!)
  constsDF$row <- constsDF$row + 1L
  # join back headers to constsDA
  constsDF <- constsDF[c("header", "row", "col", "group")]  |>
    merge(constsDA[, !names(constsDA) %in% c("group", "n"), drop = FALSE],
          all.x = T, sort = F, by = c("row", "col"))
  # we are not interested in row differences before datastart
  constsDF$row <- replace(constsDF$row,
                         constsDF$row <= headerrows[length(headerrows)],
                         headerrows[length(headerrows)] + 1L)

  # with many rows, we need a fast way to combine 2 header rows and
  # match headers to cols. Also should allow for missing headers...
  headerDF <- allCells[allCells$row == headerrows[1L] |
                         allCells$row == headerrows[2L],
                       c("character", "row", "col")] |>
    reshape(direction = "wide", timevar = "row", idvar = "col",
            v.names = "character")
  #change category for old firmware (pre 2?):
  headerDF[headerDF[[2L]] == "Sys", 2] <- "SysObs"
  headerDF$header <- paste(headerDF[[2L]], headerDF[[3L]], sep = ".") |>
    rename_header("Li6800")

  #Include datastart, this prevents that cols with no data are dropped
  #but make them NA
  dataDF <- allCells[allCells$row %in% c(datastart,datarows), ]
  dataDF$character[dataDF$row == datastart] <- NA
  # merge the double header
  idx <- match(dataDF$col, headerDF$col)
  dataDF$header <- headerDF$header[idx]
  dataDF$group <- headerDF[[2L]][idx]

  # convenient is the first row, but lets take one with formulas
  dataRow1 <- dataDF[dataDF$row == dataDF$row[!is.na(dataDF$formula)][1], ]

  # extract the formulas. I'm going to assume formulas don't change with
  # rows, except when new constants are listed between rows, but since we
  # made nice columns out of those, all formulas in a col should be
  # identical and so we use only the first row of data & consts
  # to construct R formulas
  if (extract_formula == TRUE) {

    constsRow1 <- constsDF[constsDF$row == constsDF$row[1L], ]
    # for name - address lookups:
    Row1 <- vec_rbind(dataRow1, constsRow1)[c("header", "address")]

    formDF <- dataRow1[!is.na(dataRow1$formula),
                       c("formula", "address", "header")]
    # constants become rows, there are no fixed refs anymore
    # replace fixed first. We need to strip $ before getXLrefs anyway
    formDF$formula <- stri_replace_all_fixed(formDF$formula,
                                             c("$", "TRUE()", "FALSE()"),
                                             c("",  "TRUE",   "FALSE"),
                                             vectorise_all = FALSE)

    # bind this to 2 columns containing a list with refs and derefs (names)
    formDF <- vec_cbind(formDF, t(getXLrefs(formDF$formula, Row1$address,
                                            Row1$header)))

    # now replace refs with derefs and xl_to_r functions
    formDF$fn <- stri_replace_all_regex(formDF$formula,
                                       c(paste0("\\b",
                                                unlist(formDF$refs),
                                                "\\b"),
                                         names(.gasanalyzerEnv$xl_to_r)),
                                       c(unlist(formDF$derefs),
                                         unname(.gasanalyzerEnv$xl_to_r)),
                                       vectorise_all = FALSE)

    # use a topo sort to make sure the mutate call will be ordered correctly:
    formDF <- formDF[order(tsort(formDF$address, formDF$refs)), ]
  } else
    formDF <- NA

  if (extract_units == TRUE) {
    # always keep this off:
    old_opt <- units_options("simplify")
    units_options("simplify" = NA)

    unitsDF <- allCells[allCells$row == headerrows[3L],
                        c("col", "character")]
    # replace silly things
    unitsDF$character <- stri_replace_all_fixed(unitsDF$character,
                                                c("\u207b", "\u00b2",
                                                  "\u00b9", "\u00b3", "hrs"),
                                                c("-", "2", "1", "3", "h"),
                                                vectorize_all = FALSE)
    names(unitsDF)[names(unitsDF) == "character"] <- "units"
    unitsDF$header <- dataRow1$header[dataRow1$col %in% unitsDF$col]
    units_options("simplify" = old_opt)
  } else
    unitsDF <- NA

  return(list(remarks = remarksDF, headers = headerDF, units = unitsDF,
              formula = formDF, consts = constsDF, data = dataDF))
}

#' Reads 6800 xlsx files and creates a tibble with gas-exchange data.
#'
#' The xlsx files stored by the 6800 contain measured and calculated values that
#' are read by this function and formatted in a large tibble for use with R.
#' Constants and metadata (such as calibration information) are also added as
#' columns.
#'
#' Note that values for many derived gas-exchange parameters are not
#' stored in the files, but are calculated by equations stored in the xlsx.
#' These values are 0 after importing, unless setting recalculate = TRUE. It is
#' also possible to calculate this parameters after importing using the
#' [recalculate()] function.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match. Recalculation can be disabled
#' for speed, and instead applied to the merged data using [recalculate()].
#'
#' @param filename an xlsx file containing 6800 gas-exchange data.
#' @param recalculate character string indicating whether or not to recalculate
#'   data using equations from the xlsx file.
#'
#' @returns A tibble with gas-exchange data in columns.
#'
#' @importFrom tidyxl xlsx_cells
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
#'   stri_split_fixed
#' @importFrom units set_units units_options
#' @importFrom tibble tibble as_tibble
#' @importFrom jsonify from_json
#' @importFrom tools file_path_sans_ext
#' @importFrom vctrs vec_split vec_rbind
#'
#' @seealso [recalculate()]
#'
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # read data:
#' li6800 <- read_6800_xlsx(paste0(exampledir, "//lowo2.xlsx"))
#' li6800_norecalc <- read_6800_xlsx(paste0(exampledir, "//lowo2.xlsx"),
#'   recalculate = FALSE)
#' li6800_norecalc$gasanalyzer.Equations <-
#'   list(read_6800_equations(paste0(exampledir, "//lowo2.xlsx")))
#'
#' all.equal(li6800, recalculate(li6800_norecalc), check.attributes = FALSE)
#'
#'
read_6800_xlsx <- function(filename, recalculate = TRUE) {

  DFs <- parse_6800_xlsx(filename, extract_formula = recalculate,
                         extract_units = TRUE)

  # convenience:
  reDF <- DFs[["remarks"]]
  dtDF <- DFs[["data"]]
  ctDF <- DFs[["consts"]]
  unDF <- DFs[["units"]]

  serialnumber <- reDF$X2[reDF$X1 == "Head s/n"]
  headcal <- from_json(reDF$X2[reDF$X1 == "Head cal"])

  if (length(headcal) > 1L)
    # i'll keep this as char list for now
    headcal <- rapply(headcal, as.character, how = "unlist")
  else {
    warning("Error reading head calibration from xlsx.")
  }

  filedate <- reDF$X2[reDF$X1 == "File opened"]
  if (length(filedate) == 0L) {
    warning("Error reading file date from ", filename, ".\n")
    filedate <- Sys.time()
  } else {
    # lack of tz info is cumbersome
    filedate <- as.POSIXct(filedate, tz = Sys.timezone(),
                           format = "%Y-%m-%d %H:%M:%S")
  }

  # some files have no data
  if (is.na(DFs["data"]))
    return(tibble())

  # first row (which importantly has still the header!!)
  dataRow1 <- dtDF[dtDF$row == dtDF$row[1L], ]

  constsRow1 <- ctDF[ctDF$row == ctDF$row[1L], ]
  # not sure this is really needed as the same info will be in measDF, however
  # those are not necessarily unique?
  #Row1 <- vec_rbind(dataRow1, constsRow1)[c("header", "address", "group")]

  # It may occasionally be useful to use the calculated formula values
  # the files are actually exported uncalculated by the instrument, but
  # for QA a comparison with a spreadsheet calc could be made
  formulaVector <- !is.na(dtDF$formula)
  # assume relevant formula content is numeric...
  # This slightly breaks for Const.UseDynamic
  dtDF$numeric[formulaVector] <- as.numeric(dtDF$content[formulaVector])

  #add metadata here:
  #FIXME: if cal date (firmware 2.1.11), then use that!
  factcal <- list(get_factory_cals(serialnumber, filedate)[,1])

  metadata <- list(SysObs.Filename = file_path_sans_ext(basename(filename)),
                   SysObs.Instrument = list("Li6800"),
                   SysConst.UserCal = list(headcal),
                   SysConst.FactCal = factcal)
  #NB: needs to match metadata above:
  grp <- c("SysObs", "SysObs", "SysConst", "SysConst")

  metaDF <- data.frame(row = ctDF$row[1L], header = names(metadata),
                       group = grp)
  metaDF$value <- metadata
  dtDF$isData <- TRUE
  # Following steps are costly if many rows
  measDF <- vec_rbind(dtDF, ctDF, metaDF)

  # none and "-" are sometimes used
  # one could debate whether this means 0 or NA, i'll go
  # with NA for now
  measDF$character <- replace(measDF$character,
                              measDF$character %in% c("none", "-"),
                              NA)

  exclude_meta <- seq(1, nrow(measDF) - nrow(metaDF))
  # somewhat slow:
  measDF$value[exclude_meta] <- blend(as.list(measDF$character[exclude_meta]),
                                      as.list(measDF$numeric[exclude_meta]),
                                      as.list(measDF$logical[exclude_meta]))

  # preserve row order because in a next step we are using fill for the
  # Consts rows interjected in between data rows
  measDF <- measDF[order(measDF$row), c("row", "col", "header", "value",
                                        "isData")]

  # not entirely happy with this, but seems faster than alternatives:
  measDF$rowix <- as.integer(as.factor(measDF[["row"]]))
  #fix names to unified names and remove funny chars:

  #this makes unique names at the same time
  tmp <- vec_split(measDF, measDF[["header"]])

  #since every col should have the same nr
  #of rows, cant we just complete them? here? (and save a fill later?)

  #initialize the output
  uniRows <- unique(measDF[["row"]])
  out <- rep_len(list(rep(NA, length(uniRows))), length(tmp$val) + 1)
  names(out) <- c("row", tmp$key)
  out$row <- uniRows
  #this vector stores names for filling down later:
  constsAndMeta <- rep(NA_character_, length(tmp$val))

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  #and allocate:
  for (i in seq_along(tmp$val)) {
    nm <- tmp$key[i]
    x <- tmp$val[[i]]
    out[[nm]][x$rowix] <- unlist(x$value, recursive = FALSE)

    if(is.na(x$isData[1])) {
      constsAndMeta[i] <- nm
      #only continue the loop for data:
      next
    }

    tmpu <- unDF$units[unDF$col == x$col[1L]]
    if (length(tmpu) != 0L)
      out[[nm]] <-
      # need to guard set_units as its easy to break
      tryCatch(set_units(as.numeric(out[[nm]]), tmpu,
                         mode = "standard"),
               error = function(e) {
                 warning("Failed to apply unit specified in data file: ",
                         tmpu, ", to col ", nm, "\n")
                 out[[nm]]
               })
  }
  #tibble, and drop "row":
  out <- as_tibble(out[-1L])

  # fill constants down
  constsAndMeta <- constsAndMeta[!is.na(constsAndMeta)]
  out[constsAndMeta] <- filldown(out[constsAndMeta])
  # only now drop rows that were related to Consts but had no obs (col 1)
  out <- out[!is.na(out[[1L]]), ]

  #Apertures get special treatment to remove unnecessary units:
  out["ChambConst.Aperture"] <- gsub("^([0-9.]*)(.*)$", "\\1",
                                  g0("ChambConst.Aperture", NA,
                                     envir = list2env(out)))
  # firmware bug: wrong units cf_a (fixed >=2.1.11)
  # FIXME: header says mmol should be umol, fix this better
  out["MchStatus.CFaCO2"] <- as.numeric(g0("MchStatus.CFaCO2", NA,
                                           envir = list2env(out)))

  if (recalculate) {
    eqs <- extract_6800_equations(DFs[["formula"]])
    out$gasanalyzer.Equations <- list(eqs)
    out <- recalculate(out)
  }
  out <- fixup_import(out) |>
    # add O2 uncorrected vals, this is a bit time consuming
    calculate_raw()

  units_options("simplify" = old_opt)
  #sort and return:
  out[sort_names(names(out))]

}

#' Reads 6800 text files and creates a tibble with gas-exchange data.
#'
#' The text files stored by the 6800 contain measured and calculated values that
#' are read by this function and formatted in a large tibble for use with R.
#' Constants and metadata (such as calibration information) are also added as
#' columns. Note that no recalculation of derived variables is performed,
#' although it is possible to so using [recalculate()] after importing the data.
#'
#' Multiple files can be loaded by calling the function with [lapply()] or
#' [purrr::map()] to merge multiple files. In this case, it is important
#' to ensure that the column names will match.
#'
#' @param filename an text file containing 6800 gas-exchange data.
#'
#' @returns A tibble with gas-exchange data in columns.
#'
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
#'   stri_split_fixed stri_detect_regex
#' @importFrom tibble tibble as_tibble
#' @importFrom units set_units units_options
#' @importFrom jsonify from_json
#' @importFrom tools file_path_sans_ext
#'
#' @seealso [recalculate()]
#'
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # read data
#' li6800 <- read_6800_xlsx(paste0(exampledir, "//lowo2.xlsx"))
#' li6800_txt <- read_6800_txt(paste0(exampledir, "//lowo2"))
#'
#' # compare all except equations. Note txt file reports some NAs as zero:
#' columns_to_check <- names(li6800)[!names(li6800) %in%
#'                                    c("gasanalyzer.Equations")]
#' all.equal(li6800[columns_to_check],
#'           li6800_txt[columns_to_check],
#'           tol = 0.01)
#'
read_6800_txt <- function(filename) {

  # there are various faster options, but readLines isn't that slow for
  # typical files. Note that some files seem nul-truncated
  rawfile <- readLines(filename, skipNul = TRUE)

  # find sections. search only first 1000 lines for header, and
  # assume header < 1000 lines
  rfl <- length(rawfile)
  header_start <- which(rawfile[1:min(rfl, 1000)] == "[Header]")
  if (length(header_start) == 0L) {
    warning("No [Header] section in first 1000 lines, ignoring",
            filename, ".\n")
    return(tibble())
  }

  data_start <- which(rawfile[header_start:min(rfl, header_start +
                                                 1000)] == "[Data]")
  if (length(data_start) == 0L) {
    warning("No [Data] section found, ignoring", filename, ".\n")
    return(tibble())
  }
  #stri instead of strsplit because we want to split in exactly 2:
  remarks <- stri_split_fixed(rawfile[(header_start + 1):(data_start - 1)],
                              pattern = "\t", n = 2, simplify = T,
                              omit_empty = TRUE)

  # get metadata
  bluestem <- remarks[remarks[,1] == "Console ver", 2]
  if (length(bluestem) == 0 && !grepl("Bluestem", bluestem, fixed = TRUE)) {
    warning(filename, " is not a valid 6800 file.\n")
    tibble()
  }
  filedate <- remarks[remarks[,1] == "File opened", 2]
  if (length(filedate) == 0L) {
    warning("Error reading file date from", filename, ".\n")
    filedate <- Sys.time()
  } else {
    # lack of tz info is cumbersome
    filedate <- as.POSIXct(filedate, tz = Sys.timezone(),
                           format = "%Y-%m-%d %H:%M:%S")
  }

  serialnumber <- remarks[remarks[,1] == "Head s/n", 2]
  headcal <- from_json(remarks[remarks[,1] == "Head cal", 2])

  if (length(serialnumber) > 0L && length(headcal) > 1L) {
    # i'll keep this as char list for now
    headcal <- rapply(headcal, as.character, how = "unlist")
  } else {
    warning("Error reading head calibration from", filename, ".\n")
    headcal <- NA
  }

  rawdata <- rawfile[(data_start + 1):rfl]
  # data start with a number
  # only for larger files, stringi is faster than grepl
  datarows <- stri_detect_regex(rawdata,"^[[:digit:]]+\t")
  # keep headers and consts for the moment:
  datamat <- stri_split_fixed(rawdata, pattern = "\t", simplify = TRUE)

  if (!any(datarows)) {
    warning("No data rows in ", filename, ".\n")
    return(tibble())
  }

  # remove cols that are completely empty (meanly to work around bugs
  # happening during file creation...)
  datamat <- datamat[ , colSums(datamat != "") != 0]
  # initial constants (did anybody think about how to distinguish consts
  # from remarks and other info?). I'm matching for two groups of characters
  # separated by :, but not exclusively numbers:
  # NB: they may always have Const: in the name, but that is not written
  # anywhere?
  consts_in <- stri_detect_regex(remarks[ , 1], "^(?![0-9]+).+\\:(?![0-9]+).+")
  remarks[consts_in, 1] <-
    stri_replace_all_fixed(remarks[consts_in, 1], ":", ".")

  #consts in data:
  cdata <- datamat[!datarows, 1:2]

  consts_dt <- stri_detect_regex(cdata[ ,1], "^(?![0-9]+).+\\:(?![0-9]+).+")
  cdata[consts_dt, 1] <-
    stri_replace_all_fixed(cdata[consts_dt, 1], ":", ".")
  #grp change for old firmware (pre 2?)
  datamat[1, datamat[1, ] == "Sys"] <- "SysObs"
  header <- paste(datamat[1,], datamat[2,], sep = ".")
  #Remove empty headers from our data matrix!
  nohead <- header == "."
  # BUG workaround: if the last column has data but no header, it could be due
  # to a missing Const.Custom (in mol m⁻² s⁻¹) column
  # FIXME: do more testing for this bug, it may be even worse...
  if (nohead[length(nohead)]) {
    if (all(header != "Consts.Custom")) {

      #find the custom col:
      cuscol <- which(header == "Const.Geometry")
      if (length(cuscol) > 0) {
        cuscol <- cuscol + 1
        lastcol <- length(nohead)
        #for rows where the last col is not empty
        #bring the spurious value to the end (i.e. rotate those vectors)
        badrows <- datamat[ ,lastcol] != ""
        datamat[badrows, (cuscol):lastcol] <- datamat[badrows,
                                                      c((cuscol+1):lastcol,
                                                        cuscol), drop = F]
        #but only replace from Const.Geometry onwards

        #rest of the col gets NA
        datamat[!badrows, lastcol] <- NA
        datamat[3, lastcol] <- "mol*m^-2*s^-1"
        header[lastcol] <- "Const.CustomBUG"
        nohead[lastcol] <- F
      }
    } else
      stop(filename, "is not a valid 6800 data file.")
  }
  datamat <- datamat[, !nohead]
  # these seem to mean usually NA, I hope it isn't too blunt:
  datamat <- replace(datamat, datamat %in% c("none", "-"), NA)
  header <- header[!nohead] |>
    rename_header("Li6800") |> make.unique()
  # and apply:
  colnames(datamat) <- header

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  # units (units without headers ignored!)
  unitsvec <- datamat[3, ] |>
    # remove silly things
    stri_replace_all_fixed(c("\u207b", "\u00b2", "\u00b9", "\u00b3", "hrs"),
                           c("-", "2", "1", "3", "h"), vectorize_all = FALSE)
  names(unitsvec) <- header

  #empty (but don't cut yet) headers:
  datamat[1:3, ] <- ""

  conheader <- c(remarks[consts_in, 1], cdata[consts_dt, 1]) |>
    rename_header("Li6800")
  convals <- c(remarks[consts_in, 2], cdata[consts_dt, 2])
  # put initial constants on row 1 (the header has been emptied already)
  # The data-interspaced consts go to the rows on which they are found.
  # (not on datarows, and only on those identified by consts_dt)
  conrow <- c(rep(1, sum(consts_in)),
             (seq_len(nrow(datamat)))[!datarows][consts_dt])
  # Reshape and combine with data in one go

  urow <- seq_len(nrow(datamat))
  headerall <-  c(header, conheader)
  uheaderall <- unique(headerall)
  uconheader <- unique(conheader)
  condat <- matrix(nrow = length(urow),
                   ncol = length(uheaderall),
                   dimnames = list(urow, uheaderall))
  condat[, header] <- datamat[ , header]
  condat[cbind(conrow, conheader)] <- convals

  #remaining empty cells are really NA:
  condat[condat == ""] <- NA
  # the filling and slicing may be faster on matrix, but my filldown
  # doesn't support that (yet)
  condat <- as_tibble(condat)
  #bug workaround part 2:
  if (length(condat[["Const.CustomBLC"]]) > 0) {
    #this happens sometimes, but i've no clue what the correct value should be!
    #condat[["Const.CustomBLC"]][condat[["Const.CustomBLC"]] == "True"] = "1"
    #condat[["Const.CustomBLC"]][condat[["Const.CustomBLC"]] == "False"] = "0"
    if (length(condat[["Const.CustomBUG"]]) > 0) {
      condat$Const.CustomBUG[condat$Const.CustomBUG == ""] <- NA
      condat$Const.CustomBLC <- blend(condat$Const.CustomBUG,
                                      condat$Const.CustomBLC)
      condat$Const.CustomBUG <- NULL
    }
  }
  condat[uconheader] <- filldown(condat[uconheader])
  condat <- condat[datarows,]
  condat <- units_convert(condat, unitsvec)

  condat["SysObs.Filename"] <- file_path_sans_ext(basename(filename))
  #note [[]] here
  condat[["SysConst.UserCal"]] <- list(headcal)
  #FIXME: if cal date (firmware 2.1.11), then use that!
  condat[["SysConst.FactCal"]] <- list(get_factory_cals(serialnumber,
                                                      filedate)[ , 1])
  condat[["gasanalyzer.Equations"]] <- list(NA)

  condat["SysObs.Instrument"] <- "Li6800"
  #Apertures get special treatment to remove unnecessary units:
  condat["ChambConst.Aperture"] <- gsub("^([0-9.]*)(.*)$", "\\1",
                                        g0("ChambConst.Aperture", NA_real_,
                                           envir = list2env(condat)))
  #firmware bug: wrong units cf_a (fixed >=2.1.11)
  #FIXME: may need a better fix, see above
  condat["MchStatus.CFaCO2"] <- as.numeric(g0("MchStatus.CFaCO2", NA_real_,
                                              envir = list2env(condat)))

  condat <- fixup_import(condat) |>
    # add O2 uncorrected vals, this is a bit time consuming
    calculate_raw()

  units_options("simplify" = old_opt)
  condat[sort_names(names(condat))]
}
