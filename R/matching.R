# It turns out that people make all kind of mistakes when matching. It would
# be nice to analyze and modify the match data, and recalculate gas-exchange
# parameters using this modified matching information.

# I find it hard to design an easy-to-use interface for this. This file is WIP.

#' Helper function to determine the type of matching used for gas-exchange data
#'
#' The li6800 instrument can use several different methods to calculate ref vs
#' sample offsets. It also provides detailed information about the matching in
#' the data files. However, it does not list which type of method was used to
#' calculate the matching offsets. This function attempts to use match status
#' information to determine the matching type that was used
#
#' @param me.match MchEvent.co2_match or MchEvent.h2o_match
#' @param me.at MchEvent.co2_at or MchEvent.h2o_at
#' @param ms.fl MchStatus.co2_fit_low or MchStatus.h2o_fit_low
#' @param ms.fh MchStatus.co2_fit_high or MchStatus.h2o_fit_high
#' @param a,b,c,d MchStatus.cf_<co2,h2o>_<a,b,c,d> polynomial coefficients
#' @param dlt delta value in mol fraction units defining the range for partial
#'   fits (typically maxrange/10)

#' @return A character vector with the match type
#'
#' @importFrom units drop_units
#' @keywords internal
# matchtype <- function(me.match, me.at, ms.fl, ms.fh, a, b, c, d, dlt) {
#
#     a <- drop_units(a)
#     oldc <- drop_units(me.at - me.match)
#     a_est <- drop_units(me.match) -
#         (b * oldc + c * oldc^2 + d * oldc^3)
#     #if bcd=0 we are out of range
#     ifelse(drop_units(ms.fh) == -1,
#             "point",
#             ifelse(b | c | d,
#                     ifelse(near(a, a_est, tol = 0.006), # tol because rounding
#                             ifelse(pmax(me.at - ms.fl, ms.fh - me.at) > dlt,
#                                     "full.point_adj",
#                                     "partial.point_adj"),
#                             ifelse(pmax(me.at - ms.fl, ms.fh - me.at) > dlt,
#                                     "full",
#                                     "partial")),
#                     "point.out_of_range"))
# }

#rematch <- function(conc, a, b, c, d) {
#    conc <- drop_units(conc)
    #[] tricks to maintain units of a in output

#    a[] <- drop_units(a) + (b * conc + c * conc^2 + d * conc^3)
#    a
    #me.match <- drop_units(a) +
    #a_est <- drop_units(me.match) -
    #    (b * oldc + c * oldc^2 + d * oldc^3)
#}

#' Extract match information from a tibble with gas-exchange data.
#'
#' The li6800 instrument can use several different methods to calculate ref vs
#' sample offsets. It also provided detailed information about the matching in
#' the data files. However, it does not list which type of method was used to
#' calculate the matching offsets.
#'
#' This function attempts to use the status information in a gas-exchange tibble
#' to determine the matching type that was used. By default it returns a tibble
#' with distict match information. Alternatively, it can add the match type
#' for every row in the data.
#
#' @param df a tibble with MchEvent and MchStatus fields
#' @param distinct A string specifying whether only rows with distinct values
#'   for range match coefficients excluding offset (cf_a) should be returned
#'   (distinct = "range", the default), or distinct values for all
#'   coefficients including the latest point patch (distinct = "all") or that
#'   all rows should be returned with matching information (any other value).
#
#' @return A tibble with the MchType.co2 and MchType.h2o, other matching
#'   information and (if distinct == F) possibly other columns
#'
#' @importFrom units set_units
#' @export
# extract_match <- function(df, distinct = "range") {
#
#     # delta values as defined by li6800
#     dlth2o <- set_units(2.5, "mmol mol-1")
#     dltco2 <- set_units(0.2, "mmol mol-1");
#     gas <- c("co2", "h2o")
#     # NOTE: I'm referring to these by _indices_:
#     needed_cols <- c(paste0("MchEvent.", gas, "_match"),
#                      paste0("MchEvent.", gas, "_at"),
#                      paste0("MchStatus.", gas, "_fit_low"),
#                      paste0("MchStatus.", gas, "_fit_high"),
#                      paste0(paste0("MchStatus.cf_", gas, "_",
#                                    rep(letters[1:4], each = 2))))
#     useful_cols <- c("Status.Tirga", "SysObs.time",
#                      "MchStatus.MatchCO2", "MchStatus.MatchH2O")
#
#     if (!all(needed_cols %in% colnames(df)))
#         stop("Extracting matches needs MchEvent.* and MchStatus.* columns.")
#
#     # SysObs.time and Tirga are facultative but allow distinct output to have
#     # a timestamp and temperature for the first measurement they were used for
#     df %>% {if (distinct == "all")
#         arrange(., across(any_of("SysObs.time"))) %>%
#             select(., any_of(useful_cols),
#                    all_of(needed_cols)) %>%
#             #distinct _match and _cf:
#             distinct(., across(all_of(needed_cols[c(1,2, 9:16)])),
#                      .keep_all = T)
#         else if (distinct == "range")
#             arrange(., across(any_of("SysObs.time"))) %>%
#             select(., any_of(useful_cols),
#                    all_of(needed_cols)) %>%
#             #distinct cf_[bcd], but filter cf=0:
#             filter(!if_all(needed_cols[11:16],~ .x == 0)) %>%
#             distinct(., across(all_of(needed_cols[11:16])),
#                      .keep_all = T)
#         else . } %>%
#         #use syms for co2 and h2o cols
#         mutate(MchType.co2 = matchtype(!!!syms(needed_cols[c(T,F)]),
#                                        dltco2),
#                MchType.h2o = matchtype(!!!syms(needed_cols[c(F,T)]),
#                                        dlth2o),
#                .before = 1)
# }

#' Modify match information in a gas-exchange data.
#'
#' It can be useful to alter match information in an existing gas-exchange
#' tibble. For example when forgetting to enable range matching or when a
#' correct range match was done after the measurement.
#'
#' It is not so easy to come up with a convenient interface for this. This
#' function allows to replace range matches in the data using two returned
#' tibbles from matchExtract: one listing the match coefficients to replace,
#' one containing the replacement coefficients. If more specific match
#' adjustments are desired, manual replacement and recalculation is necessary.
#'
#' After successful replacement rematch is called to recalculate
#' concentrations and all gas-exchange parameters
#' TODO: this last thing is not done yet!
#'
#' @param df a gas-exchange tibble with MchEvent and MchStatus fields
#' @param src a tibble with indexes to match the rows in df that are to
#'   be replaced by dst
#' @param dst a tibble with MchStatus coefficients to replace src. Must be
#'   the same number of rows as src, or just one row.
#' @param adj point adjusted range matches are used by default. If dst contains
#'   match information without point adjustment, this can be set to FALSE.
#
#' @return A recalculated gas-exchange tibble with replaced match information.
#
#' @importFrom units set_units
#' @importFrom stringi stri_subset_regex
#' @export
# modify_match <- function(df, src, dst, adj = T) {
#     gas <- c("co2", "h2o")
#     # NOTE: I'm referring to these by _indices_:
#     needed_cols <- c(paste0("MchEvent.", gas, "_match"),
#                      paste0("MchEvent.", gas, "_at"),
#                      paste0("MchStatus.", gas, "_fit_low"),
#                      paste0("MchStatus.", gas, "_fit_high"),
#                      paste0(paste0("MchStatus.cf_", gas, "_",
#                                    rep(letters[1:4], each=2))))
#
#     if (!all(needed_cols %in% colnames(df)))
#         stop("Input dataframe requires MchEvent.* and MchStatus.* columns.")
#
#     if(nrow(src) != nrow(dst) & nrow(dst) != 1)
#         stop(paste("If dst has more than one row, src and dst need to have ",
#                    "the same number of rows."))
#
#     cf_co2 <- stri_subset_regex(c(colnames(dst)),
#                                 "^MchStatus.cf_co2_[bcd]$")
#     cf_h2o <- stri_subset_regex(c(colnames(dst)),
#                                 "^MchStatus.cf_h2o_[bcd]$")
#
#     co2.ops <- (length(cf_co2) == 3) & any(length(cf_h2o) == c(0, 3))
#     h2o.ops <- (length(cf_h2o) == 3) & any(length(cf_co2) == c(0, 3))
#
#     if (!(co2.ops | h2o.ops))
#         stop(paste("Coefficients b, c and d for either CO2 or H2O",
#                    "must be present in dst."))
#
#     #TODO: if MchType fails we could call matchtype()
#     if (adj == F) {
#         if ((co2.ops & !all(c("MchStatus.cf_co2_a", "MchType.co2") %in%
#                             colnames(dst))) |
#             (h2o.ops & !all(c("MchStatus.cf_h2o_a", "MchType.h2o") %in%
#                             colnames(dst))))
#             stop(paste("If point adjustment is not used, dst for CO2 or ",
#                  "H2O must contain all 4 coefficients and a MchType column."))
#         if (co2.ops) cf_co2 <- c(cf_co2, "MchStatus.cf_co2_a")
#         if (h2o.ops) cf_h2o <- c(cf_h2o, "MchStatus.cf_h2o_a")
#     }
#
#     #NOTE: cf_* contains dupes, but this apparently works
#     dst <- dst[grepl(colnames(dst), "(cf_co2)|(cf_h2o)")]
#
#     replace.cols <- colnames(dst)
#     src <- cbind(src, rename_with(dst, ~ paste0(.x,".dst")))
#
#     # left_join our replacement table to the df so we match
#     # on all columns in src. Then replace cols in dst using
#     # coalesce when !NA in *.dst. Map over all dst cols and
#     # finally remove *.dst cols.
#     # not sure this is the best way to do it
#     df %>%
#         left_join(src) %>%
#         mutate(as_tibble(Map(function(x, y) coalesce(y, x),
#                              .[, replace.cols],
#                              select(., paste0(replace.cols, ".dst"))))) %>%
#         select(-ends_with(".dst"))
#
# }
