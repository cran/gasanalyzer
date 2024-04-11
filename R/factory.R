# factory calibrations

#' Import instrument-specific factory calibration files from a folder.
#'
#' The factory calibration of the 6800 can be used to calculate
#' concentrations from raw values. The calibration files are found
#' in subfolders of /home/licor/.factory on the instrument. They can be copied
#' to a computer and imported to the package configuration files using this
#' method.
#'
#' The function will also load the calibration into the package environment,
#' where they can be retrieved by [get_factory_cals()].
#'
#' @details
#' This method assumes the files are named with serial number and calibration
#' date, separated by an underscore.
#'
#' @param folder A folder where calibration files are to be found.
#' @param keep Copies valid calibration files to a package-specific
#'   configuration folder. Will result in automatic import of the data in the
#'   future. Will overwrite files.
#' @returns Calibration data is stored in the package environment.
#'
#' @importFrom jsonify from_json
#' @importFrom utils type.convert modifyList
#' @importFrom stringi stri_split_fixed stri_replace_all_regex
#' @importFrom stats setNames
#' @importFrom tools R_user_dir
#' @seealso [get_factory_cals()]
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#'
#' # show calibration data
#' get_factory_cals()
#'
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # show calibration data
#' get_factory_cals()
#'
import_factory_cals <- function(folder = tools::R_user_dir("gasanalyzer",
                                                         which = "config"),
                              keep = FALSE) {

  defaultDir <- eval(formals(import_factory_cals)$folder)

  if (!dir.exists(defaultDir))
    dir.create(defaultDir, recursive=TRUE)

  if (!dir.exists(folder))
    stop(folder, " not found.")

  #assumes files are named with serial_date:
  files <- list.files(folder, pattern = "*._.*")

  # this defines what names and types need to be present in a factory
  # calibration file so that we can safely rely on it in the calculations
  # add more here when needed.
  irn <- list(a1 = NA, a2 = NA, a3 = NA, a4 = NA, a5 = NA, z = NA, s1 = NA,
              s2 = NA, s3 = NA)
  bblist <- list(ch = "1.4", cx = "0.9", hx0 = "0.96", hx1 = "0.00190476",
                 pcorr_c1 = NA, pcorr_c2 = NA,
                 pcorr_c3 = NA, pcorr_h1 = NA,
                 pcorr_h2 = NA, pcorr_h3 = NA)
  example <- unlist(list(serialnumber = NA, caldate = "01 Jan 1970",
                         cal = list(time = "0", type = "Factory"),
                         bb = bblist,
                         irga_a = list(Xhc = "0.001", Xch = "-0.0001",
                                       co2 = irn, h2o = irn[-c(4,5)]),
                         irga_b = list(Xhc = "0.001", Xch = "-0.0001",
                                       co2 = irn, h2o = irn[-c(4,5)])),
                    recursive = T)
  nme <- names(example)

  sn <- stringi::stri_replace_all_regex(files, "^([\\w-]+)(_[\\d-\\w]+)", "$1")
  sn2 <- c(unique(sn), "default")
  out <- setNames(vector("list", length(sn2)), sn2)
  out[["default"]][["0"]] <- example

  if (length(files) > 0) {
    for (i in seq_along(files)) {
      y <- from_json(file.path(folder, files[i]))
      if (length(y) != 0) {
        # validate json names and classes:
        flaty <- unlist(y)
        ms <- nme[!(nme %in% names(flaty))]
        if (length(ms) > 0) {
          warning(files[i], " is missing: ", paste(ms, collapse = ", "),
                  ".\n Calibration file invalid.")
          # we want to prevent keeping such files:
          files[i] <- NA
        } else {
          # note we merge identical sn and distinguish by the timestamp
          out[[sn[i]]][[y[["cal"]][["time"]]]] <- flaty[nme]
        }

      }
    }

    files <- files[!is.na(files)]
    if (length(files) > 0) {
      # Ordering using filename would only have allowed ordering on the date,
      # now we can use the actual timestamp:
      # make sure that files are in the correct order:
      out <- lapply(out[lengths(out) !=0],
                    function(x) {x[order(as.numeric(names(x)),
                                         decreasing = T)]})
      if (keep & !missing(folder)) {
        message("Copying ", length(files), " calibration files to ",
                defaultDir, ".")
        file.copy(file.path(folder, files), defaultDir, overwrite = T)
      } else
        packageStartupMessage("Loaded ", length(files), " factory ",
                              "calibrations from ", folder, ".")
    }
  }
  # assign for later use
  assign("factory", out, .gasanalyzerEnv)

}

#' Returns a matrix with factory calibration information for given instrument
#' serial numbers and calibration dates.
#'
#' The factory calibration of the 6800 can be used to calculate
#' concentrations from raw values. If calibration information is available in
#' the package environment it can be retrieved by this method.
#'
#' @param sn a character vector with an instrument serial number. If named, the
#'   names are kept in the output.
#' @param datetime a POSIXct time vector indicating the latest possible time for
#'   the calibration data that is to be returned. If no calibration before
#'   datetime is found, the oldest available calibration is returned.
#' @returns A character matrix with factory calibration data. If no datetime is
#'   provided, the newest calibration is returned.
#'
#' @details The datetime option can be used to make sure that newer calibration
#'   files are not used in combination with older datafiles.
#'
#' @importFrom jsonify from_json
#' @importFrom utils type.convert
#' @importFrom stats setNames
#'
#' @export
#' @examples
#' exampledir <- system.file("extdata//", package = "gasanalyzer")
#'
#' # import factory cals for example data:
#' import_factory_cals(exampledir)
#'
#' # show calibration data for a specific instrument serial numbers, closest to
#' # the current time:
#' get_factory_cals(sn = "68H-422400", datetime=Sys.time())
get_factory_cals <- function(sn = NULL, datetime = NULL) {

    cals <- get0("factory", .gasanalyzerEnv)
    idx <- match(sn, names(cals), nomatch = 0)
    if (!is.null(sn)) cals <- setNames(cals[idx], sn[idx != 0])
    if (length(cals) < 1 ) return(matrix(NA))

    # Both mapply and sapply simply the output to a matrix, but only if possible
    # This should nevertheless be save as unlist is called first....
    if (!is.null(datetime)) {
      # every entry in datetime must have a matching entry in the list
      mapply(function(x,y) { tdiff <- as.numeric(y) - as.numeric(names(x))
                             # if we fail to find a date before the measurement
                             # simply return the earliest (length(x) calibration
                             unlist(x[c(which.min(tdiff[tdiff>0]),
                                        length(x))][[1]])
                           }, cals, datetime)
    } else
      sapply(cals, function(x) unlist(x[[1]]))
}

# I leave this here as a WIP in case cal sheets get fixed.

#vecex <- function(v1, v2) {
#  idx <- which(v1 == v2[1L])
#  v2l <- length(v2) - 1L
#  idx[sapply(idx, function(i) all(v1[i:(i + v2l)] == v2))]
#}

#cal_from_pdf <- function(filename) {
  # download?
#if (!requireNamespace("pdftools", quietly=TRUE))
# stop("Parsing pdf calibration sheets requires the pdftools package.")

#  pdfdata <- pdftools::pdf_data(filename)[[1]]
  # did any body bother to look at the dump of machine generated data
  # in a PDF that is mangled to unparsable xml?

  # select the small font
#  pdfdata <- pdfdata[pdfdata$height < 5,]
  # assume we have 3 columns (lets hope that is reliable)
#  pdfdata$groups <- kmeans(pdfdata$x, centers=3)$cluster
#  pdfdata <- pdfdata %>%
#      group_by(groups) %>%
#      arrange(y, .by_group = T) %>%
#      pull(text)
  # indexed because they might actually occur several times in this weird doc...
#  datastart <- last(vecex(pdfdata, c("<licor>","<li6850>","<factory>")))
#  dataend <- vecex(pdfdata, c("</factory>","</li6850>","</licor>"))[1]
#  pdfdata <- pdfdata[datastart:dataend]

#}

