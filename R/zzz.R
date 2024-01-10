utils::globalVariables(c("."))

#' Globals and initializations
#' @import units
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @noRd
.onLoad <- function(libname, pkgname) {

  package_ns <- parent.env(environment())
  #options we use
  options("gasanalyzer.calibration.warning" = TRUE)
  units::units_options("negative_power" = TRUE, simplify = NA)
  units::remove_unit("steps")
  units::install_unit("steps")
  #permille:
  units::remove_unit("\U2030")
  units::install_unit(symbol = "\U2030", def = "0.001", name = "permille")
  #empty units look better than 1 for plotting
  units::remove_unit("\U200B")
  units::install_unit(symbol = "\U200B", def = "1", name = "dimensionless")

  # NOTE: units has some options that appear to change its behavior
  # which will lead to different (wrong) results. Simplify seems to be
  # the main offender

  assign("unity", units::set_units(1e6, "ppm", mode = "standard"),
         envir = package_ns)

  assign("gasconstant", units::set_units(8.314462618, "J*K^-1*mol^-1",
                                  mode = "standard"), envir = package_ns)

  load_internal_data()
  import_factory_cals()
}
