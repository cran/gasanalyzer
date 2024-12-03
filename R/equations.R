#' Create a list of equations for recalculating gasanalyzer data.
#'
#' This function creates a list of equations that can be used to recalculate
#' gas-exchange data by passing the resulting object to the [recalculate()]
#' method. Various \code{useflags} can be defined to tune the equations.
#' In addition, custom equations can be defined as arguments. Note that
#' the calculations may fail if commons are missing in the gas-exchange data.
#'
#' @param useflags character vector with the type of equations to create
#'   (such as c("li6800", "gfs3000")). Leave empty to obtain the default set.
#'   An unknown flag returns an empty list, and a warning listing all valid
#'   flags.
#' @param ... custom equations. the arguments must tagged function expressions.
#'   Note that the function body must be wrapped in curly brackets. The tags
#'   will be matched against the names of a data frame when applying the return
#'   value with [recalculate()].
#'
#' @returns A list of language objects with equations
#'
#' @importFrom utils getParseData modifyList
#' @importFrom stats setNames ave
#' @importFrom stringi stri_split_fixed stri_list2matrix stri_wrap
#' @importFrom vctrs vec_duplicate_detect vec_duplicate_id
#' @seealso [read_6800_equations()]
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#'
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # read data from a txt file:
#' li6800 <- read_6800_txt(file.path(exampledir, "lowo2"))
#'
#' # passing an invalid flags shows which flags are valid:
#' \donttest{create_equations("help")}
#'
#' # create a default set of gas-exchange equations, for the 6800, but overwrite
#' # the default calculation of leaf light absorption with a custom value:
#' Eqs <- create_equations(c("default", "li6800"), LeafQ.alpha = \() {0.86})
#'
#' #apply:
#' li6800_recalc <- recalculate(li6800, Eqs)
#'
#' li6800$LeafQ.alpha
#' li6800_recalc$LeafQ.alpha
create_equations <- function(useflags = "default", ...) {

  if(!is.null(useflags) &&
     (typeof(useflags) != "character" ||
      any(make.names(useflags, unique = TRUE) != useflags)))
    stop("`useflags` must be a character value or vector of syntactically ",
         "valid unique names.")

  ueqs <- unlist(.gasanalyzerEnv$vars$fn[!is.na(.gasanalyzerEnv$vars$fn)])
  # splitted names:
  neqs <- stringi::stri_split_fixed(names(ueqs), ".", simplify = FALSE)
  validflags <- NULL
  fnRowIdx <- vector("logical", length = length(neqs))

  #first loop: get valid flags, add deps to useflags and get fn rows
  for (i in seq_along(neqs)) {
    x <- neqs[[i]]
    nn <- length(x)
    validflags <- union(validflags, x[3:(nn - 1)])

    if (x[[nn]] == "deps" && ueqs[[i]] != "") {
      if (all(x[3:(nn-1)] %in% useflags)) {
        useflags <- union(useflags, ueqs[[i]])
      }
    } else if (x[[nn]] == "fn" && length(ueqs[[i]]) != 0)
      fnRowIdx[i] <- TRUE

  }

  notaflagIdx <- !(useflags %in% validflags)
  if (any(notaflagIdx)) {
    wtext <- paste0(validflags[!(validflags %in% c(""))],
                           collapse = ", ") |>
      stri_wrap(initial = "Valid flags are: ") |> paste(collapse = "\n")
    warning("\n", wtext, "\n")
    useflags <- useflags[!notaflagIdx]
  }

  #second loop get equations
  for (i in which(fnRowIdx)) {
    x <- neqs[[i]]
    if (!all(x[3:(length(x) - 1)] %in% useflags)) {
      fnRowIdx[i] <- FALSE
    }
  }

  if (!any(fnRowIdx == TRUE) || length(useflags) == 0) {
    warning("No valid equations found!\n")
    return(list(NULL))
  }

  ueqs <- ueqs[fnRowIdx]
  nms <- stri_list2matrix(neqs[fnRowIdx], byrow = TRUE)
  names(ueqs) <- paste(nms[ ,1], nms[ ,2], sep = ".")
  #dups need removal. Assign a score based on how many flags match. The default
  #flag does not count towards the score, so its never the preferred option:
  scr <- nms[ , 3:(ncol(nms))] %in% useflags[!(useflags %in% "default")] |>
    as.numeric() |> matrix(nrow = nrow(nms)) |> rowSums()
  #find max score per duplicated item:
  dup.id <- vctrs::vec_duplicate_id(nms[ , 1:2])
  non.dup <- as.logical(ave(scr, dup.id, FUN = \(x) x == max(x)))
  if (length(unique(dup.id)) < length(which(non.dup))) {
    dups <- unique(nms[vec_duplicate_detect(nms[non.dup, 1:2]), 3])
    stop("Cannot simultaneously apply: ", paste(dups, collapse = " and "), ".")
  }
  ueqs <- ueqs[non.dup]

  tokens <- lapply(ueqs,
                   function(x) {
                     y <- as.list(x)[[1]]
                     attributes(y) <- NULL
                     y
                   })

  custom.eqs <- as.list(substitute(...()))
  #or: rlang::enexprs(..., .ignore_empty = "all", .named = T,
  #                           .ignore_null = "all", .homonyms = "last")
  #this would allow !!, so might be better
  # loop to filter out non-funcs
  if(any(make.names(names(custom.eqs), unique = TRUE) != names(custom.eqs)))
    stop("`...` must use syntactically valid and unique names.")

  custom.eqs <- lapply(custom.eqs, function(x) {
    if (is.function(eval(x)) && length(x) == 4) x[[3]] else NULL
  })

  if (!all(lengths(custom.eqs) > 1)) {
    warning("Non-function arguments specified in ... will be ignored.\n")
  }

  # this also conveniently removes NULL values (made by !is.function above)
  tokens <- modifyList(tokens, custom.eqs)

  # we now loop again to get dependencies of all equations:
  refs <- vector("list", length = length(tokens))
  for (i in seq_along(tokens)) {
    pd <- getParseData(parse(text = tokens[i], keep.source = TRUE))
    sbls <- unique(pd$text[pd$token == "SYMBOL"])
    refs[[i]] <- sbls[sbls %in% names(tokens)]
  }

  tokens <- tokens[order(tsort(names(tokens), refs))]
  #rlang::new_quosures(as_quosures(tokens, parent.env(environment())))
  tokens[["gasanalyzer.UseFlags"]] <- useflags
  tokens
}

#' Modify an existing list of equations with specific user-specified equations.
#'
#' This method allows replacing a specific equations in a list with
#' custom versions. Although it is possible to add custom equations using
#' [create_equations()], it can be useful to modify existing sets. It can also
#' be used to modify equations imported from an `xlsx` file.
#'
#' @param eqs a list of calls for recomputing `gasanalyzer` equations.
#' @param ... custom equations. the arguments must tagged function expressions.
#'   The tags will be matched against the equation list specified in eqs, and
#'   matching expressions will be replaced. Additional expressions will be
#'   added to the list. Note that the function body must be wrapped in curly
#'   brackets.
#'
#' @returns A modified list of calls containing equations to recalculate
#'   `gasanalyzer` data.
#'
#' @importFrom utils getParseData modifyList
#' @importFrom stats setNames ave
#' @importFrom stringi stri_split_fixed stri_list2matrix stri_wrap
#' @importFrom vctrs vec_duplicate_detect vec_duplicate_id
#' @seealso [read_6800_equations()]
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#'
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # read data from a txt file:
#' li6800 <- read_6800_txt(file.path(exampledir, "lowo2"))
#'
#' # create a default set of gas-exchange equations, for the Li-6800:
#' Eqs <- create_equations(c("default", "li6800"))
#'
#' # replace the value for the leaf light absorptance:
#' Eqs <- modify_equations(Eqs, LeafQ.alpha = \() {0.86})
#'
#' # apply:
#' li6800_recalc <- recalculate(li6800, Eqs)
#'
#' li6800$LeafQ.alpha
#' li6800_recalc$LeafQ.alpha
modify_equations <- function(eqs, ...) {

  if(is.null(eqs) || !inherits(eqs, "list") ||
      length(eqs) < 1 || !is.call(eqs[[1]]))
    stop("`eqs` must be a list with calls for recomputing gasanalyzer data.")

  #FIXME: code from here essentially duplicates create_equations, should be
  #refactored
  custom.eqs <- as.list(substitute(...()))

  if(any(make.names(names(custom.eqs), unique = TRUE) != names(custom.eqs)))
    stop("`...` must use syntactically valid and unique names.")

  custom.eqs <- lapply(custom.eqs, function(x) {
    if (is.function(eval(x)) && length(x) == 4) x[[3]] else NULL
  })

  if (!all(lengths(custom.eqs) > 1)) {
    warning("Non-function arguments specified in ... will be ignored.\n")
  }

  # this also conveniently removes NULL values (made by !is.function above)
  eqs <- modifyList(eqs, custom.eqs)

  # we now loop again to get dependencies of all equations:
  refs <- vector("list", length = length(eqs))
  for (i in seq_along(eqs)) {
    pd <- getParseData(parse(text = eqs[i], keep.source = TRUE))
    sbls <- unique(pd$text[pd$token == "SYMBOL"])
    refs[[i]] <- sbls[sbls %in% names(eqs)]
  }

  eqs <- eqs[order(tsort(names(eqs), refs))]
  #rlang::new_quosures(as_quosures(eqs, parent.env(environment())))
  eqs

}

#' Recalculate gas-exchange data based on a set of equations.
#'
#' The recalculation uses equations in a list of quosures provided as argument.
#' This list can be obtained from [create_equations()] or
#' [read_6800_equations()].
#
#' @param df A data frame or an extension thereof (e.g. a tibble).
#' @param eqs a list of quosures that define how the df will be altered.
#'
#' @returns A tibble with recalculated columns as specified by the eqs
#'   argument
#'
#' @importFrom vctrs vec_split vec_rbind
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
#' exampledir <- system.file("extdata", package = "gasanalyzer")
#' # import factory calibration for example data:
#' import_factory_cals(exampledir)
#'
#' # read data:
#' li6800 <- read_6800_xlsx(file.path(exampledir, "lowo2.xlsx"))
#'
#' # recalculate using xlsx equations:
#' li6800 <- recalculate(li6800)
#'
#' # recalculate using gasanalyzer default equations for the li6800:
#' li6800_ge <- recalculate(li6800, create_equations(c("default", "li6800")))
#'
#' # the difference is that units have been enforced using gasanalyzer, which
#' # has been recorded in a column:
#' all.equal(li6800, li6800_ge[names(li6800)], tol = 0.01)
recalculate <- function(df, eqs = NULL) {

  if (length(nrow(df)) < 1 || nrow(df) < 1L || ncol(df) < 2) {
    warning("Input df is too smal, returning as-is.\n")
    return(df)
  }

  if (length(eqs) == 0 || !inherits(eqs, "list") ) {
    if (length(df[["gasanalyzer.Equations"]]) == 0) {
      warning("No valid equations provided or found in df. Returning as-is.\n")
      return(df)
    }
    df <- lapply(vec_split(df, df[["gasanalyzer.Equations"]])[[2]],
           function(x) { .recalculate(x) }) |>
      do.call(vec_rbind, args = _) |>
      #the rbind may change the order
      (\(x) { x[order(x[["SysObs.Obs"]]), ]})()
  } else {
    df <- .recalculate(df, eqs)
  }
  as_tibble(df[sort_names(names(df))])
}

#' Internal function handling the recalculations
#'
#' @importFrom units drop_units
#' @noRd
.recalculate <- function(df, eqs = NULL) {

  if (length(eqs) == 0) {
    #df should have been splitted on equations, so take first element here:
    eqs <- df[["gasanalyzer.Equations"]][[1]]
    if (length(eqs) == 0 || !is.list(eqs)) {
      warning("No equations provided or found for all rows. Returning as-is.\n")
      return(df)
    }
  }

  # always keep this off:
  old_opt <- units_options("simplify")
  units_options("simplify" = NA)

  package_ns <- parent.env(environment())

  useflags <- eval(eqs[["gasanalyzer.UseFlags"]])
  if (length(useflags) == 0) useflags <- "default"

  df <- fixup_import(df, useflags, FALSE)
  eqs[["gasanalyzer.UseFlags"]] <- NULL
  # get names _after_ setting UseFlags to NULL
  nms <- names(eqs)
  ununit <- isFALSE(as.logical(eval(eqs[["gasanalyzer.UseEqUnits"]])))
  if (ununit) {
    #add cols if not there to prevent setting a NULL class below.
    #we sadly can't know units for these...
    df[nms[!(nms %in% names(df))]] <- NA_real_
    ununit.df <- drop_units(df)

    for (i in seq_along(eqs)) {
      tmp <- df[[nms[[i]]]]
      try(tmp <- eval(eqs[[i]], ununit.df, package_ns))
      ununit.df[[nms[[i]]]] <- tmp
      class(tmp) <- class(df[[nms[[i]]]])
      attr(tmp, "units") <- attr(df[[nms[[i]]]], "units", T)
      df[[nms[[i]]]] <- tmp
    }
  } else {
    for (i in seq_along(eqs)) {
      # evaluate within the package environ to use @ and others
      try(df[[nms[[i]]]] <- eval(eqs[[i]], df, package_ns))
    }
  }
  units_options("simplify" = old_opt)
  df
}

