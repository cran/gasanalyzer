#' Short-hand expression to handle units.
#'
#' Since it is cumbersome to write long equations that need manual tweaking to
#' avoid unit errors, this provides a simple way to enforce or drop units.
#'
#' @param x an object to apply units to.
#' @param y a unit to convert x to. . (a period) will drop units.
#'
#' @returns an object with units metadata (or with units removed if y==".").
#'
#' Note that it is possible to assign units using a character variable, but be
#' careful when such variables redefine existing unit names!!
#'
#' @importFrom methods slot
#' @importFrom units set_units
#' @noRd
"@" <- function(x, y = ".") {
    ys <- substitute(y)
    if (!isS4(x)) {
      if (exists(deparse(ys), parent.frame(1)) && is.character(y)) ys <- y[1]
      if (ys == ".") as.numeric(x) else {
        storage.mode(x) <- "numeric"
        set_units(x, ys, mode="standard")
      }
    } else
      slot(x, ys)
}

#' Find first non-missing element.
#'
#' like [dplyr::coalesce], except that vctrs::vec_equal_na got updated to
#' ignore NA in lists; So blend instead replaces NA also in lists.
#'
#' @param ... One or more vectors. Vectors are recycled to a common length
#'   and cast to a common type.
#' @returns A vector with the same type and size as in ..., with missing
#'   elements filled.
#'
#' @importFrom vctrs vec_size vec_cast_common vec_recycle_common vec_slice
#' @noRd
blend <- function(...) {
  args <- list(...)

  n_args <- vec_size(args)
  if (n_args == 0L) return(NULL)

  args <- do.call(vec_cast_common, args)
  args <- do.call(vec_recycle_common, args)

  out <- args[[1L]]
  args <- args[-1L]

  for (arg in args) {
    is_na <- is.na(out)
    if (!any(is_na)) break
    vctrs::vec_slice(out, is_na) <- vctrs::vec_slice(arg, is_na)
  }
  out
}

#' Fill missing values with previous value.
#'
#' Like dplyr::fill, but for downwards direction, based on ideas described here:
#' stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value
#'
#' @param df A data frame or tibble.
#'
#' @returns A data frame or tibble where NA values have been substituted with
#'   the latest non-NA value.
#' @noRd
filldown <- function(df) {
  for (i in seq_along(df)) {
    x <- df[[i]]
    ind <- which(!is.na(x))
    if(is.na(x[1])) ind <- c(1, ind)
    df[[i]] <- rep(x[ind], times = diff(c(ind, length(x) + 1)))
  }
  df
}

#' Apply fixes and sanitize imported data.
#'
#' Imported gas-exchange data may miss various columns, or have incorrect data
#' types or units. This function will take a dataframe and apply types and units
#' according to the pre-specified rules. It will also add missing columns with
#' default values.
#'
#' For internal use.
#'
#' @param df A data frame or tibble
#' @param useflags character vector with flags that indicate which columns get
#'  added/checked based on an internal list
#' @param units boolean to indicate whether units should be applied to all
#'  columns or only to the default values
#' @returns A data frame or tibble with the corrections applied
#'
#' @importFrom units set_units
#' @noRd
fixup_import <- function(out, useflags = "default", units = TRUE) {

  useflags <- useflags[useflags %in% .gasanalyzerEnv$vars$useflag]
  vars <- .gasanalyzerEnv$vars[.gasanalyzerEnv$vars$useflag %in% useflags, ]
  # first set defaults (this makes them character but that will be fixed later):
  defaultsNms <- vars$fullname[!is.na(vars$default)]
  defaultsVals <- vars$default[!is.na(vars$default)]
  #NB: this adds missing cols, but not replaces NA!
  out[defaultsNms] <- with(out, mget(defaultsNms, ifnotfound = defaultsVals))
  out[["SysObs.Obs"]] <- with(out, get0("SysObs.Obs",
                                        ifnotfound = seq.int(1, nrow(out))))

  #units only for defaults to speed up things
  if (units == FALSE)
    nms <- defaultsNms
  else
    nms <- colnames(out)

  #prevent units being converted:
  type.convert.units <- function(x, ...) {x}

  #applying units is not vectorized:(
  for (x in nms) {
    un <- c(vars$units[vars$fullname == x], NA)[1]
    if (length(un) > -10) {
      out[[x]] <- switch(un,
                         "posix"    = as.POSIXct(as.numeric(out[[x]]),
                                                 tz = Sys.timezone(),
                                                 origin = "1970-01-01"),
                         "numeric" = ,
                         "character" = ,
                         "logical" = ,
                         "double" = ,
                         "integer" = {
                           #storage.mode preserves units!
                           storage.mode(out[[x]]) <- un
                           out[[x]]
                         },
                         `NA` = type.convert(out[[x]], as.is = TRUE,
                                             trylogical = FALSE),
                         "list" = ,
                         "ignore" = out[[x]],
                         { storage.mode(out[[x]]) <- "numeric"
                           set_units(out[[x]], un, mode = "standard")
                         }
      )
    } else
        out[[x]] <- type.convert(out[[x]], as.is = TRUE, trylogical = FALSE)
  }

  out
}

#' Orders gas-exchange names in a standard way.
#'
#' Imported gas-exchange data may not be in a unified order. A large
#' number of variables and constants are known to the package and this allows
#' to sort the columns in a consistant way. Unknown columns are sorted
#' according to the group they are in.
#'
#' @param nms A character vector of names.
#' @returns A character vector with sorted names.
#'
#' @noRd
sort_names <- function(nms) {
  #this wastes some microseconds because we have group names already
  #but in theory, they could have been altered by rename_header:
  grps <- stringi::stri_split_fixed(nms, ".", simplify = TRUE)[ , 1]
  #presort on group:
  idx <- grps %in% .gasanalyzerEnv$vars$group
  tosort <- nms[idx]
  nms[idx] <- tosort[order(factor(grps[idx], ordered = TRUE,
                                  levels = unique(.gasanalyzerEnv$vars$group)))]
  #now sort on name:
  idx <- nms %in% .gasanalyzerEnv$vars$fullname
  tosort <- nms[idx]
  nms[idx] <- tosort[order(factor(nms[idx], ordered = TRUE,
                                  levels = .gasanalyzerEnv$vars$fullname))]

  nms
}

#' Rename a header vector.
#'
#' Translate the vendor names to universal and typeable names
#' The lookup table is stored in a specific environment. The names are not
#' necessarily unique.
#'
#' @param x vector with names.
#' @param nm a character value with the instrument name.
#'
#' @returns a renamed vector with names.
#'
#' @importFrom stringi stri_replace_all_fixed stri_detect_fixed
#' @noRd
rename_header <- function(x, nm) {

  idx <- x %in% .gasanalyzerEnv$vars[[nm]]
  x[idx] <- .gasanalyzerEnv$vars[["fullname"]][match(x,
                                                     .gasanalyzerEnv$vars[[nm]],
                                                     0)]
  # only replace funny symbols in the rest:
  x[!idx] <-   stri_replace_all_fixed(x[!idx],
                                      .gasanalyzerEnv$sympairs[c(T, F)],
                                      .gasanalyzerEnv$sympairs[c(F, T)],
                                      vectorise_all = F)

  idx2 <- stri_detect_fixed(x[!idx], ".")
  x[!idx][!idx2] <- paste(nm, x[!idx][!idx2], sep = ".")
  x
}


#' Expand a data frame with all possible combinations of the values in a column.
#'
#' For sensitivity analyses, it is useful to permutate the values in a single
#' column, whilst keeping all other values constant. After creating such a
#' permutation, [recalculate()] should be used to analyze the effect of the
#' change in the column of interest. If the effect of changes in multiple
#' columns is to be analyzed, this function can be called in series.
#'
#' @param df a dataframe or tibble
#' @param ... a name-value pair. The name gives the name of the column in the
#'   input that is to be changed. The value is a vector specifying all values
#'   that are desired in the output. For every value in this vector, all other
#'   rows are duplicated.
#'
#' @returns a data frame containing all possible combinations of the input df
#'   and the vector specified in ...
#'
#' Note that the units and classes of the columns in the input data frame are
#' applied to the replacement values. Unexpected behavior may occur when
#' providing incompatible classes or units.
#'
#' @importFrom vctrs vec_rep
#' @importFrom stringi stri_paste
#' @export
#'
#' @examples
#' example <- system.file("extdata", "6400-testfile", package = "gasanalyzer")
#'
#' # read data:
#' li6400 <- read_6400_txt(example)
#'
#' # expand the data frame for a range of leaf areas, and recalculate the data:
#' li6400 <- permutate(li6400, Const.S = seq(1, 8)) |>
#'   recalculate(create_equations(c("default", "li6400")))
#'
#' if (interactive()) {
#'   require(units)
#'   require(graphics)
#'
#'   # observe that changing the leaf area enclosed in the chamber would have a
#'   # nonlinear effect on the rate of photosynthesis:
#'   aggregate(list(A = li6400$GasEx.A), list(Area = (li6400$Const.S)), mean) |>
#'     plot()
#' }
permutate <- function(df, ...) {

  pCol <- as.list(substitute(...()))

  #in theory, we could do more at the same time....
  if (length(pCol) > 1) {
    warning("Please specify a single column to permutate.")
    pCol <- pCol[[1]]
  } else if (length(pCol) == 0) {
    warning("No column to permutate specificied.")
    return(df)
  }

  pColNm <- names(pCol)
  if (length(pColNm) == 0 || !pColNm %in% names(df))
    stop("Column names in ... must already be defined in df.")

  n <- nrow(df)
  un <- get_units(df[[pColNm]])
  #would be nice if uc accepts vectors as well:
  #to support lists this is convoluted:
  tmp <- eval(pCol[[1]])
  edf <- data.frame(b = rep(NA, length(tmp)))
  edf[["a"]] <- tmp
  newVals <- units_convert(edf, c(a = un))[["a"]]

  df <- vec_rep(df, length(newVals))
  #NB: for named lists, this gives the rows names that may be identical
  #this can be useful, but also may be confusing?
  df[[pColNm]] <- rep(newVals, each = n)
  df[["gasanalyzer.PermutateLabel"]] <- rep(names(tmp), each = n)
  #add a label for what was permutated, not sure: list/mat/char?
  #note stri_paste allows avoiding trailing separators
  df[["gasanalyzer.Permutate"]] <- stri_paste(pColNm,
                                               df[["gasanalyzer.Permutate"]],
                                               sep=",", ignore_null = TRUE)

  df


}

#' Converts Comment fields to oxygen levels
#'
#' GFS software adds comments when you change the oxygen setting in the menu
#' But it does not add an O2 column to the data. We use the comment vector
#' to generate an O2 vector. If there's no comment referring to O2 on the
#' first line, we assume it is 21%.
#'
#' It does not yet fill rows between comments.
#'
#' @param Comment a vector with comments.
#' @returns a vector with O2 concentrations.
#'
#' @importFrom stringi stri_match_last
#' @noRd
comment_to_oxygen <- function(Comment)
{
    re <- "(?<O2>(?:\\d{1,3}\\.\\d{0,5})(?=% O2)|(?:\\d{1,3})(?=% O2)|(?:air))"
    O2Col <- stri_match_last(Comment,
                             regex = paste0("(Calculate H2O-data with )",
                                            re))[,"O2"]
    O2Col <- replace(O2Col, O2Col == "air", "21")
    # assume no comment means 21% O2
    if (is.na(O2Col[1])) O2Col[1] <- "21"
    as.numeric(O2Col)
}


#' Deparse units, or date or time into a string.
#'
#' Typically used to generate headers saving data.
#'
#' @param data a tibble or df to process.
#' @return a vector with the units as strings.
#'
#' @importFrom units deparse_unit
#' @noRd
get_units<-function(data) {
    switch(class(data)[1],
           "units"    = deparse_unit(data),
           "POSIXct"  = "posix",
           class(data)[1]
    )
}

#' Piecewise interpolation.
#'
#' The approx method in stats takes all x and y values into account.
#' Instead, it can be usefull to do a piecewise interpolation between
#' every two subsequent points.
#'
#' @param x,y vectors of x and y values.
#' @param z a logical vector which values in y need to be interpolated.
#' @returns a modified vector y. Elements for which z == TRUE have been
#'   interpolated.
#'
#' @importFrom stats approx
#' @noRd
pwapprox <- function(x, y , z) {

    known_i <- which(!z)
    n <-  length(known_i)-1

    for(i in 1:n) {
        xvals <- x[c(known_i[i], known_i[i+1])]
        yvals <- y[c(known_i[i], known_i[i+1])]

        xout <- x[known_i[i]:(known_i[i+1]-1)]
        y[known_i[i]:(known_i[i+1]-1)] <- approx(xvals, yvals, xout,rule = 2)$y
    }
    y
}

#' Convert character columns to specified units or other types as appropriate.
#'
#' Note that unitsvec is a named vector with the names corresponding to the
#' column names of df. Empty or NA elements get passed on to type.convert.
#'
#' @param df a dataframe to convert.
#' @param unitsvec a named vector with units that are passed on to set_units.
#'
#' @returns a df with columns using the appropriate units or storage types.
#'
#' @importFrom units set_units
#' @noRd
units_convert <- function(df, unitsvec) {

  unitsvec <- unitsvec[unitsvec != ""]
  for (x in names(df)) {
    df[[x]] <- switch(unitsvec[x],
           "posix"    = as.POSIXct(as.numeric(df[[x]]),
                                   tz = Sys.timezone(),
                                   origin = "1970-01-01"),
           "numeric" = ,
           "character" = ,
           "logical" = ,
           "double" = ,
           "integer" = {
             storage.mode(df[[x]]) <- unitsvec[x]
             df[[x]]
            },
           `NA` = df[[x]],
           "list" = df[[x]],
           {
              storage.mode(df[[x]]) <- "numeric"
              tryCatch(set_units(df[[x]], unitsvec[x], mode = "standard"),
                       error = function(e) {
                         warning("Failed to apply unit ",
                                 unitsvec[x], " to col ", x, "\n")
                         df[[x]]
                       })
            }
           )
  }

  df
}

#' An implementation of a topological sort.
#'
#' An implementation of Kahn's algorithm. It is used here to order formulas so
#' they are executed in an order that makes sure prerequisite values are before
#' values depending on them.
#
#' @param nm a character vector with items.
#' @param deps a list of vectors, indicating dependencies for every nm.
#'
#' @return a numeric vector with elements indicating the order of nm in deps
#' Later elements depend on earlier elements. The output can be used with
#' [base::order()] or [dplyr::arrange()] to rearrange nm.
#'
#' @importFrom stats setNames
#'
#' @noRd
tsort <- function(nm, deps) {
  refs <- union(as.vector(unlist(deps)), nm)
  deps <- setNames(deps, nm)
  S <- setdiff(refs, nm) #entries that depend on nothing
  k <- TRUE
  while(k) {
    k <- FALSE
    for(x in setdiff(nm, S)) {
      R <- c(S, x)
      if(length(setdiff(deps[[x]], R)) == 0) {
        S <- R
        k <- TRUE
      }
    }
  }

  if(length(S) < length(refs))
    warning(sprintf(", circular dependencies found: \n%s",
                    paste("", setdiff(refs, S), sep = "", collapse = "\n")))

  # output order of nm in s
  match(nm, S)

}

#' Allows calling a variable that may not exist.
#'
#' @param x a variable.
#' @param ifnotfound the return value when x does not exists. Defaults to 0.
#' @param un units of the return variable.
#' @param na_replace treat NA values as not existing. Defaults to FALSE.
#'
#' @noRd
g0 <- function(x, ifnotfound = 0, un = NULL, na_replace = FALSE,
               envir = parent.frame(1)) {

  x <- get0(substitute(x), envir = envir, ifnotfound = ifnotfound)

  if (na_replace) x[is.na(x)] <- ifnotfound
  if (length(un) == 0) return(x)
  #this doesn't change units class objects, so better than as.numeric:
  storage.mode(x) <- "numeric"
  set_units(x, un, mode = "standard")
}

#' Helper function to rename columns
#'
#' @param df a data.frame or tibble
#' @param old_names character vector with names of columns that are potentially
#'   in df, and should be renamed
#' @param new_names the new names for the columns. Needs to be a character
#'   vector with the same length as old_names
#'
#' @returns a data.frame or tibble with changed names
#' @noRd
rename_cols <- function(df, old_names, new_names) {
  nms <- colnames(df)
  ex_nms <- old_names[old_names %in% nms]

  colnames(df)[match(ex_nms, nms)] <- new_names[match(ex_nms, old_names)]
  df
}

#' Use the water mol fractions to back-calculate O2 concentrations for
#' the GFS-3000.
#'
#' The calculation is based on a table provided by the company.
#'
#' @param H2O reference H2O mol fraction.
#' @param wa sample H2O mol fraction.
#' @param mp uncorrected difference between reference and sample.
#' @param zp GFS zero-point correction (match-offset).
#'
#' @returns the O2 correction factor.
#'
#' @importFrom units set_units drop_units
#' @noRd
gfs_calc_o2 <- function(H2O, wa, mp, zp) {
  # The provided table seems to use a piecewise linear fit of
  # several H2O levels. at 500 ppm the values are identical to 1000 ppm H2O.
  # Similarly, for 20% O2, values from 20.942 seem to be reused.
  # I used a linear fit (fit is near perfect) of the factor against O2 to find
  # the slopes at key H2O levels
  o2cal <- data.frame(
    h2o = c(0, 1000, 5000, 10000, 15000, 20000, 30000, 45000, 75000),
    o = c(1, 0.99, 0.966, 0.955, 0.949, 0.946, 0.941, 0.933, 0.928),
    s = c(0, 4.775096e-04, 1.623532e-03, 2.148791e-03, 2.435297e-03,
          2.57855e-03,  2.817305e-03, 3.199313e-03, 3.438067e-03))
  H2O <- as.numeric(set_units(H2O, "umol/mol"))
  H2Of <- as.numeric(set_units(wa, "umol/mol")) /
    (as.numeric(mp) - as.numeric(zp) + H2O)
  if (any(H2Of > 1.01 | H2Of < 0.90, na.rm = T))
    warning ("Unlikely values for calculating O2 concentration.")

  below <- findInterval(H2O, o2cal$h2o, all.inside = T)
  above <- below + 1
  x1 <- (o2cal$h2o[above] - H2O)
  x2 <- (H2O - o2cal$h2o[below])
  p1 <- o2cal$o[below] * x1 + o2cal$o[above] * x2
  p2 <- o2cal$s[below] * x1 + o2cal$s[above] * x2

  round((H2Of * (o2cal$h2o[above] - o2cal$h2o[below]) - p1) / p2)

}

#' Calculates a GFS-3000 specific correction factor for the water vapor delta
#' based on the background water and oxygen levels.
#'
#' The factor is based on a table provided by the company. Note that
#' no data is available for O2 levels above 21% (factor assumed 1).
#'
#' @param H2O reference H2O mol fraction.
#' @param O2 reference O2 mol fraction.
#' @returns the O2 correction factor.
#'
#' @importFrom units set_units drop_units
#' @noRd
gfs_o2_factor <- function(H2O, O2) {
  # The provided table seems to use a piecewise linear fit of
  # several H2O levels. at 500 ppm the values are identical to 1000 ppm H2O.
  # Similarly, for 20% O2, values from 20.942 seem to be reused.
  # I used a linear fit (fit is near perfect) of the factor against O2 to find
  # the slopes at key H2O levels
  o2cal <- data.frame(
    h2o = c(0, 1000, 5000, 10000, 15000, 20000, 30000, 45000, 75000),
    o = c(1, 0.99, 0.966, 0.955, 0.949, 0.946, 0.941, 0.933, 0.928),
    s = c(0, 4.775096e-04, 1.623532e-03, 2.148791e-03, 2.435297e-03,
          2.57855e-03,  2.817305e-03, 3.199313e-03, 3.438067e-03))

  O2 <- as.numeric(set_units(O2, "%"))
  H2O <- as.numeric(set_units(H2O, "umol/mol"))

  if (any(O2 > 21, na.rm = T))
    warning ("No calibration data available for O2 levels above 21%.")

  below <- findInterval(H2O, o2cal$h2o, all.inside = T)
  above <- below + 1
  x1 <- (o2cal$h2o[above] - H2O)
  x2 <- (H2O - o2cal$h2o[below])
  p1 <- o2cal$o[below] * x1 + o2cal$o[above] * x2
  p2 <- o2cal$s[below] * x1 + o2cal$s[above] * x2

  ((O2 * p2 + p1) / (o2cal$h2o[above] - o2cal$h2o[below])) |>
    replace((H2O < 0) | (O2 > 21), 1)

}

#' Modified Arden Buck equation
#'
#' This is currently only used for solving the EB for temperature and gbw
#' simultaneously.
#'
#' @param Temp Temperature as double representing K
#' @param Pres Pressure as double representing kPa
#' @returns Saturation vapor pressure in mol fraction
#'
#' @noRd
sat_vapor_pressure <- function(Tp, Ps) {
  Tp <- Tp - 273.15
  EF <- 1 + 1e-3 * (0.72 + Ps * (0.0320 + 5.9e-6 * Tp^2))
  (EF * (0.61121 * exp((18.678 - Tp / 234.5) * (Tp / (257.14 + Tp))))) / Ps

}

#' This function can be passed to nleqslv to solve the energy balance for
#' both temperature and boundary layer conductance. It uses the method
#' described by Bristow 1987. The Buck equation is used to calculate
#' the saturation vapor pressure.
#'
#' In the future, it could be extend to also provide a jacobian, which
#' would vastly speed up the method.
#'
#' @param Tleaf a vector with the initial estimate of the sample temperature,
#'  in units
#' @param Tair a vector with the air temperature, in units
#' @param Twall a vector with the chamber wall temperature, in units
#' @param Pcham a vector with the chamber pressure, in units
#' @param H2Os a vector with the water vapor concentration in the chamber, in
#'   mol fraction units
#' @param E a vector with the measured transpiration rate in the chamber, in
#'   units
#' @param RHi a vector with the relative humidity inside the chamber, in units
#'   or numeric
#' @param Rabs a vector with the absorbed short wave radiation of the sample,
#'   in units
#' @param lambda a vector with the latent heat of evaporation of water, in
#'   J per mol units
#' @param eps a vector with the emmisivity of the sample, in units or numeric
#' @param ra_rv a vector with the ratio between heat resistance and H2O
#'   diffusion resistance of the air, in units or numeric
#' @param Cpm a vector with the molar heat capacity of the air in the chamber,
#'   in J per mol per K units
#'
#' @returns a vector with the estimated sample temperature in Kelvin units
#'
#' @importFrom units set_units
#' @importFrom nleqslv nleqslv
#' @noRd
energy_balance <- function(Tleaf, Tair, Twall, Pcham, H2Os, E, RHi, Rabs,
                           lambda, asH, eps, ra_rv, Cpm) {
  # we dispense with units for speed, but be careful to keep all
  # EB components in W/m2
  Tl <- as.numeric(set_units(Tleaf, "K"))
  Ta <- as.numeric(set_units(Tair, "K"))
  Tw4 <- as.numeric(set_units(Twall, "K"))^4
  Pa <- as.numeric(set_units(Pcham, "kPa"))
  wa <- as.numeric(set_units(H2Os, "1"))
  E <- as.numeric(set_units(E, "mol*m^-2*s^-1"))
  l <- as.numeric(set_units(lambda, "J*mol^-1"))
  es <- as.numeric(set_units(eps, "1")) * 5.67e-08
  cc <- as.numeric(set_units(ra_rv, "1"))
  a <- as.numeric(set_units(asH, "1"))
  cp <- as.numeric(set_units(Cpm, "J*mol^-1*K^-1"))
  Ra <- as.numeric(set_units(Rabs, "W*m^-2"))
  RHi <- as.numeric(set_units(RHi, "1"))

  model <- function(x) {
    Tf <- x[c(T, F)]
    gb <- x[c(F, T)]
    wi <- sat_vapor_pressure(Tf, Pa) * RHi
    # planar leaf: a = 2 (convection and radiation on both sides)
    # but Ra only from one side. E also happens from both sides
    # hypo leaf: a=2, but E would be from one side? However, this is not
    # called for leaves with stomata
    # round leaf: a=1, meaning the gb here is a total gbw
    EB <- Ra + a * es * (Tw4 - Tf^4) + a * cp * (Ta - Tf) * gb / cc - l * E
    Ev <- a * gb * (wi - wa) / (1 - (wi + wa) / 2) - E
    c(EB, Ev)
  }
  # pass Tair and 4 mol/m2/s as initial conditions
  ss <- nleqslv(c(rbind(Tl, 4)), model, control = list(ftol = 1e-4))
  data.frame(gb = set_units(ss$x[c(F, T)], "mol*m^-2*s^-1"),
             Tf = set_units(ss$x[c(T, F)], "K"))

}


#' Render gasanalyzer variables or values using mathematical notation.
#'
#' The argument is converted to a plotmath expression that can be used using
#' `plot` or `ggplot2`. If there is no known plotmath expression defined,
#' the argument is returned as is.
#'
#' @param varname a character list or vector argument with variable names.
#' @param use_units whether or not to append default units to the resulting.
#'   expression. It is better to rely on the units in the actual data, which is
#'   handled automatically by newer ggplot2 versions.
#' @param val a value to display rather than the variable name. Needs to be of
#'   the same length as varname and if not a character vector it should be
#'   coercible to one. NAs are replaced with the plotmath
#'   expressions.
#' @param ... options passed on to make_unit_label.
#'
#' @returns a list of plotmath expressions.
#'
#' @importFrom units make_unit_label as_units
#' @seealso plotmath
#' @export
#'
#' @examples
#' # make labels
#' lbls <- var2label(c("GasEx.Ci", "GasEx.A"), use_units = TRUE)
#' print(lbls)
#'
#' # plot
#' plot(1, type = "n", xlab = lbls[[1]], ylab = lbls[[2]])
#' # add temperature as title, removing [] from the units:
#' title(main = var2label("Meas.Tleaf",use_units = TRUE, val = 25,
#'   group = c("", ""))[[1]])
#'
var2label <- function(varname, use_units = FALSE, val = NULL, ...) {
  nn <- length(varname)
  if (length(val) == 0)
    val <- rep(NA, nn)
  else if (length(val) != nn)
    stop("varname and val need to be of equal length, or val needs to be NULL.")

  if (length(use_units) == 1)
    use_units <- rep(use_units, nn)
  else if (length(use_units) != nn)
    stop("use_units needs to be of length 1 or of the length as varname.")

  #cols that match names
  mtch <- match(varname, .gasanalyzerEnv$vars$fullname)
  #symbols (might be NA!), blend in vals
  symbs <- blend(as.character(val), .gasanalyzerEnv$vars$plotmath[mtch])
  idx <- which(!is.na(mtch))
  varname <- as.list(varname)

  un <- .gasanalyzerEnv$vars$units[mtch]
  # export of ud_is_parseable would be nice (+vectorized):
  un[un %in% c("numeric", "string", "list", "logical",
               "integer", "posix") | is.na(un)] <- list("\U200B")
  varname[idx] <- lapply(idx, \(x) {
    if (use_units[x])
      make_unit_label(symbs[[x]], as_units(un[[x]]), ...) |> as.expression()
    else
      str2expression(symbs[[x]])
  })

  varname
}
