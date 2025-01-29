#' Create a list of equations for recalculating gasanalyzer data.
#'
#' This function creates a list of equations that can be used to recalculate
#' gas-exchange data by passing the resulting object to the [recalculate()]
#' method. Various `useflags` can be defined to tune the equations.
#' In addition, custom equations can be defined as arguments. Note that
#' the calculations may fail if commons are missing in the gas-exchange data.
#'
#' @param useflags character vector with the type of equations to create
#'   (see Details). Leave empty to obtain the default set. An unknown flag
#'   returns an empty list, and a warning listing all valid flags.
#' @param ... custom equations. the arguments must tagged function expressions.
#'   Note that the function body must be wrapped in curly brackets. The tags
#'   will be matched against the names of a data frame when applying the return
#'   value with [recalculate()].
#'
#' @returns A list of language objects with equations
#'
#' @details
#' The `useflags` argument currently supports several pre-defined sets
#' of equations.
#' ## Gas-Exchange and fluorescence models
#' \describe{
#'   \item{`default`}{
#'     A set with the most commonly-used derived quantities in the gas-exchange
#'     (`GasEx`) and chlorophyll fluorescence (`FLR`) categories. The
#'     equations are described in detail by von Caemmerer and Farquhar (1981)
#'     and LI-COR Biosciences Inc (2022).
#'   }
#'   \item{`GoffGratch1946`}{
#'     Replaces or adds equations related to the calculation of
#'     the saturated water pressure of the leaf (`GasEx.SVPleaf`) and chamber
#'     air (`GasEx.SVPcham`). The calculation is based on that described by
#'     Goff and Gratch (1946) and only takes the temperature into account.
#'   }
#'   \item{`Buck1981`}{
#'     Replaces or adds equations for (`GasEx.SVPleaf`) and (`GasEx.SVPcham`).
#'     based on the description by Buck (1981). Takes temperature and
#'     pressure into account.
#'   }
#'   \item{`Buck1996`}{
#'    Replaces or adds equations for (`GasEx.SVPleaf`) and (`GasEx.SVPcham`)
#'    based on the description by Buck (1996). Takes temperature and
#'    pressure into account.
#'   }
#'   \item{`cuticular_conductance`}{
#'     Replaces the equations related to CO2 and H2O conductance
#'     and substomatal CO2 (`GasEx.gtw`; `GasEx.gsw`; `GasEx.gtc`; `GasEx.Ci`)
#'     with versions that take into account cuticular conductance
#'     (Márquez et al., 2021, 2023). Requires manually specifying the
#'     cuticular conductance to water (`Const.gcw`) and CO2 (`Const.gcc`).
#'   }
#'   \item{`boundary_conducance`}{
#'     Replaces or adds equations related to the boundary
#'     layer conductance (`GasEx.gbw`) and leaf temperature derived from
#'     the energy balance equations (`GasEx.TleafEB`) with a version that
#'     is valid for a sample with no stomata. This is typically used when
#'     estimating `GasEx.gbw` using filter paper (Parkinson, 1985). Note
#'     that the leaf thermocouple (`Meas.Tleaf`) is used to estimate air
#'     temperature (and should therefore not touch the sample). This can be
#'     overridden by a custom equation. The sample emissivity (`Const.eps`),
#'     the ratio between heat and water conductance for the chamber air
#'     (`Const.ra_rv`) and the fraction of the area exchanging radiative and
#'     sensible heat (`Const.asH`) can be adjusted. For solving the energy
#'     balance equations, the method described by Bristow (1987) is used.
#'   }
#'   \item{`gm_fluorescence`}{
#'     Adds derived variables for mesophyll conductance
#'     (`FLR.gm`) and chloroplast CO2 mole fractions (`FLR.Cc`) based
#'     on gas exchange and chlorophyll fluorescence (Harley et al., 1992). It is
#'     strongly recommended to first calibrate the electron transport rates
#'     reported in the `FLR.ETR` column. In addition, the equations require
#'     adding columns for the respiration rate in the light (`Const.RL`),
#'     and the CO2 photo-compensation point (`Const.GammaStar`).
#'   }
#' }
#' ## Isotope models
#' \describe{
#'   \item{`d13C`}{
#'     Adds derived variables related to the stable carbon isotope
#'     discrimination model for C3 plants (Farquhar and Cernusak, 2012;
#'     Evans and von Caemmerer, 2013) aimed at estimating the mesophyll
#'     conductance `d13C.gm` and chloroplast CO2 mole fraction
#'     (`d13C.Cc`). Requires additional columns with data on the carbon
#'     isotope composition in sample and reference air
#'     (`d13CMeas.delta13CO2s`; `d13CMeas.delta13CO2r`) and in the air
#'     where the plants were grown (`d13CConst.delta13CO2g`), and values
#'     for respiration in the light (`Const.RL`) and the CO2
#'    photo-compensation point (`Const.GammaStar`).
#'   }
#'   \item{`d13C_dis`}{
#'     Requires the `d13C` flag, but modifies the
#'     modeled carbon isotope discrimination `d13C.Deltai` and
#'     `d13C.gm` using the assumption that the carbon pools of respiration
#'     and assimilation are disconnected (as described by Busch et al., 2020).
#'   }
#'   \item{`d13C_e_Busch2020`}{
#'     Requires the `d13C` flag, but modifies the calculation
#'     of the effective respiratory fractionation (`d13C.ep`) to
#'     better take into account the effect of the growth conditions
#'     (Busch et al. 2020). Additionally requires a value for the observed
#'     discrimination against 13CO2 under growth conditions
#'     (`d13CConst.Deltag`).
#'   }
#' }
#'
#' ## Instrument related calculations
#' These implement calculations as close as possible to those used in the
#' instrument firmware and are typically related to the specific configurations
#' or instrument design. In addition, some methods are provided for
#' recalculating low-level instrument variables.
#' \describe{
#'   \item{`li6400`}{
#'     Definitions for boundary layer conductance, temperature, and light
#'     related variables (`GasEx.gbw`; `GasEx.Rabs`;
#'     `GasEx.TairCnd`; `LeafQ.alpha`) specific to  LI-6400 / LI-6400XT
#'     instruments (LI-COR Biosciences Inc, 2011).
#'   }
#'   \item{`li6800`}{
#'     Definitions for boundary layer conductance, light, leakage,
#'     temperature and light related variables (`GasEx.gbw`;
#'     `GasEx.Rabs`; `GasEx.TairCnd`; `LeafQ.alpha`;
#'     `LeafQ.Qin`; `LeafQ.Conv`; `Leak.Fan`;
#'     `Leak.CorrFact`; `GasEx.Ca`) specific to the LI-6800 instrument
#'     (LI-COR Biosciences Inc, 2022).
#'   }
#'   \item{`ciras4`}{
#'     Light absorptance (`LeafQ.alpha`) and the conversion between
#'     photon flux density and irradiance (`LeafQ.Conv`) is calculated by
#'     taking into account the effect of the different light sources that can be
#'     used with the CIRAS-4 instrument (PP Systems, 2024).
#'   }
#'   \item{`gfs3000`}{
#'     Boundary layer conductance (`GasEx.gbw`) and light sensor
#'     (`LeafQ.Qin`) calculations adjusted for the default chamber of the
#'     GFS-3000 instrument (Heinz Walz GmbH, 2019).
#'   }
#'   \item{`gfs3000_light_bot`}{
#'     Requires gfs3000, but modifies `LeafQ.Qin` to indicate that the
#'     bottom light sensor of the default chamber was used to quantify the
#'     light intensity incident on the leaf (Heinz Walz GmbH, 2019).
#'   }
#'   \item{`match`}{
#'     Takes a previously stored offset between sample and reference
#'     analyzers into account when recalculating water and CO2 mole fractions
#'     (`Meas.H2Os`; `Meas.CO2s`) as described by LI-COR Biosciences
#'     Inc (2022), Heinz Walz GmbH (2019) and PP Systems (2024). The corrected
#'     mole fractions are already stored in the data and therefore this
#'     calculation is usually not required. However, these equations are needed
#'     when recalculating mole fractions from lower level data (see the
#'    `raw` and `O2_correction` flags).
#'   }
#'   \item{`raw`}{
#'     Recalculates CO2 and H2O mole fractions from such low-level variables
#'     (`Meas.CO2r`; `Meas.H2Or`; `Meas.CO2a`; `Meas.H2Oa`)
#'     as described by LI-COR Biosciences Inc (2022). Currently only
#'     implemented for the LI-6800 because low-level instrument data are
#'     required. Requires storing raw data and the availability of factory
#'     calibration files. Requires and enables the `match` set.
#'   }
#'   \item{`O2_correction`}{
#'     Recalculates CO2 and H2O mole fractions (`Meas.CO2r`;
#'     `Meas.H2Or`;  `Meas.CO2a`;  `Meas.H2Oa`) at a potentially
#'     different oxygen concentration ( `Const.Oxygen`). Currently only
#'     implemented for the LI-6800 (LI-COR Biosciences Inc, 2022) and GFS-3000
#'     (K. Siebke, Heinz Walz GmbH, personal communication). For the LI-6800,
#'     this requires loading of factory calibration files
#'     ([import_factory_cals()]). Requires and automatically enables the
#'     `match` flag.
#'   }
#' }
#' @references
#' * Bristow KL. 1987. On solving the surface energy balance equation for
#' surface temperature. Agricultural and Forest Meteorology **39**:49–54.
#' * Buck AL. 1981. New equations for computing vapor pressure and enhancement
#' factor. Journal of Applied Meteorology **20**:1527–1532.
#' * Buck AL. 1996. Buck research CR-1A user’s manual. Boulder, CO: Buck
#' Research Instruments.
#' * Busch FA, Holloway-Phillips M, Stuart-Williams H, Farquhar GD. 2020.
#' Revisiting carbon isotope discrimination in C3 plants shows respiration
#' rules when photosynthesis is low. Nature Plants **6**:245–258.
#' <https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf>
#' * von Caemmerer S, Farquhar GD. 1981. Some relationships between the
#' biochemistry of photosynthesis and the gas exchange of leaves.
#' Planta **153**:376–387.
#' * Evans JR, von Caemmerer S. 2013. Temperature response of carbon isotope
#' discrimination and mesophyll conductance in tobacco. Plant, Cell &
#' Environment **36**:745–756.
#' * Farquhar GD, Cernusak LA. 2012. Ternary effects on the gas exchange of
#' isotopologues of carbon dioxide. Plant, Cell & Environment
#' **35**:1221–1231.
#' * Goff JA, Gratch S. 1946. Low-pressure properties of water from -160° F to
#' 212° F. Transactions of the American Society of Heating and Ventilating
#' Engineers **52**:95–121.
#' * Harley PC, Loreto F, Di Marco G, Sharkey TD. 1992. Theoretical
#' considerations when estimating the mesophyll conductance to CO2 flux by
#' analysis of the response of photosynthesis to CO2. Plant Physiology
#' **98**:1429–1436.
#' * Heinz Walz GmbH. 2019. Portable gas exchange fluorescence system GFS-3000.
#' Handbook of operation. 9th edition.
#' <https://www.walz.com/files/downloads/manuals/gfs-3000/GFS-3000_Manual_9.pdf>
#' * LI-COR Biosciences Inc. 2011. Using the LI-6400/LI-6400XT portable
#' photosynthesis system. Version 6.2.
#' <https://www.licor.com/env/support/LI-6400/manuals.html>
#' * LI-COR Biosciences Inc. 2022. Using the LI-6800 portable photosynthesis
#' system. Version 2.1.
#' <https://www.licor.com/env/support/LI-6800/manuals.html>
#' * Márquez DA, Stuart-Williams H, Farquhar GD. 2021. An improved theory for
#' calculating leaf gas exchange more precisely accounting for small fluxes.
#' Nature Plants **7**:317–326.
#' * Márquez DA, Stuart-Williams H, Cernusak LA, Farquhar GD. 2023. Assessing the
#' CO2 concentration at the surface of photosynthetic mesophyll cells.
#' New Phytologist **238**:1446–1460.
#' * Parkinson KJ. 1985. A simple method for determining the boundary layer
#' resistance in leaf cuvettes. Plant, Cell & Environment **8**: 223–226.
#' * PP Systems. 2024. CIRAS-4 portable photosynthesis system. Operation manual.
#' Amesbury, MA: PP Systems. Version 1.3.
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

