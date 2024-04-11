#' Simple convenience function for translating XL LEFT
#
#' @param txt a character vector
#' @param nr number of characters from start of txt to keep
#' @return a character vector with some characters removed
#'
#' @importFrom stringi stri_sub
#' @noRd
left <- function(txt, nr = 1L) { stri_sub(txt, to = max(nr, 0)) }

#' Simple convenience function for translating XL RIGHT
#
#' @param txt a character vector
#' @param nr number of characters from start of txt to keep
#' @return a character vector with some characters removed
#'
#' @importFrom stringi stri_sub
#' @noRd
right <- function(txt, nr = 1L) { stri_sub(txt, from = -max(nr, 0)) }

#' Simple convenience function for translating XL radians
#
#' @param deg a numeric vector with degrees
#'
#' @return a numeric vector with radians
#'
#' @importFrom units set_units
#' @noRd
radians <- function(deg)
    as.numeric(set_units(set_units(deg, "degrees"), "radians"))

#' Simple convenience function for translating XL degrees
#
#' @param rad a numeric vector with radians

#' @return a numeric vector with degrees
#'
#' @importFrom units set_units
#' @noRd
degrees <- function(rad)
    as.numeric(set_units(set_units(rad, "radians"), "degrees"))

#' Simple convenience function for translating XL ceiling
#
#' @param x a numeric vector with numbers to round up
#' @param accuracy accuracy of number to round up to
#' @return a numeric vector with rounded numbers
#'
#' @noRd
ceiling_xl <- function(x, accuracy) {
    ceiling(x / accuracy) * accuracy
}

#' Simple convenience function for translating XL PI()
#
#' @return pi
#'
#' @noRd
pi_xl <- function() return(pi)

#' Helper function to get cell addresses out of XL equations.
#'
#' Translate addresses to variable names used in R.
#' Vectorizes the tidyxl::xlex call.
#'
#' @param frm a character vector with XL formulas
#' @param deps a character vector with addresses
#' @param nm a character vector with variable names
#'
#' @return a list with 2 elements: a vector with all addresses and a vector with
#' names that they match to.
#'
#' @importFrom tidyxl xlex
#'
#' @noRd
getXLrefs <- Vectorize(function(formula, addr, nm) {
    t <- xlex(formula)
    rf <- unique(t$token[t$type == "ref"])
    list(refs = rf, derefs = nm[match(rf, addr)])
}, vectorize.args="formula", USE.NAMES = F)

