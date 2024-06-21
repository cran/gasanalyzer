#' Gas-exchange equations.
#'
#' Currently this is only used internally to generate package provided equations
#' during initialization. I might expand it in the future to provide more
#' help on equations.
#'
#' Perhaps add descriptions: `r calcs()$Leak.Fan$default$desc`
#'
#' Or equations:
#' \deqn{\mathrm{Leak.Fan} = \max(0, \frac{\mathrm{LeakConst.fan\_a} +
#'   \mathrm{LeakConst.fan\_b} \cdot \mathrm{Meas.Fan\_speed}}
#'   {1 + \mathrm{LeakConst.fan\_c} \cdot \mathrm{Mean.Fan\_speed}}
#'   \frac{\mathrm{Meas.Pa}\cdot
#'   \mathrm{LeakConst.fan\_d}}{\mathrm{Meas.Tair}+273.15}) }
#'
#' \deqn{\mathrm{Leak.CorrFact} = \begin{cases}
#'   1 &\text{if } \mathrm{Leak.Leak} \cdot \mathrm{LeakConst.leak\_wt}
#'   >= \mathrm{Leak.Fan} \\
#'   \frac{\mathrm{Leak.Fan}}{\mathrm{Leak.Fan} - \mathrm{Leak.Leak} \cdot
#'   \mathrm{LeakConst.leak\_wt}} &\text{if } \mathrm{Leak.Leak} \cdot
#'   \mathrm{LeakConst.leak\_wt}
#'   < \mathrm{Leak.Fan} \end{cases}}
#'
#' @returns A list of supported computations.
#'
#' @noRd
calcs <- function() {

  # using "with" because functions in the list are for non-standard eval.
  df <- data.frame()

  out <- with(df, list(
    Leak.Fan = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(Leak.Fan, NA_real_, "\U00B5mol*s^-1")
        }
      ),
      li6800 = list(
        desc = "Li6800 specific equation.",
        deps = "",
        fn = function() {
          # licor txt deviates slightly from licor xlsx
          pmax(
            0,
            (LeakConst.CFaFan + LeakConst.CFbFan *
               Meas.FanSpeed@rpm@.) / (1 + LeakConst.CFcFan *
                                          Meas.FanSpeed@rpm@.) *
              Meas.Pa@"kPa"@. / Meas.Tair@"K"@. *
              LeakConst.CFdFan
          )@"\U00B5mol*s^-1"
        }
      )
    ),
    Leak.CorrFact = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { g0(Leak.CorrFact, 1, NULL) }
      ),
      li6800 = list(
        desc = "Li6800 specific version.",
        deps = "",
        fn = function() {
          ifelse(Leak.Leak * LeakConst.LeakWt >= Leak.Fan, 1,
                 (Leak.Fan / (
                   Leak.Fan - Leak.Leak * LeakConst.LeakWt))@.)
        }
      )
    ),
    Leak.CorrFactPct = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          (Leak.CorrFact - 1) * 100@"%"
        }
      )
    ),
    LeafQ.Qin = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(LeafQ.Qin, Meas.QambIn, "\U00B5mol*m^-2*s^-1")
        }
      ),
      li6800 = list(
        desc = "Li6800 specific version.",
        deps = "",
        fn = function() {
          QConst.fQambIn * g0(Meas.QambIn, 0, "\U00B5mol*m^-2*s^-1") +
            QConst.fQambOut * g0(Meas.QambOut, 0, "\U00B5mol*m^-2*s^-1") +
            QConst.fQflr * g0(FlrLS.Q, 0, "\U00B5mol*m^-2*s^-1") *
              (1@"1" - g0(FlrLS.ffarred, 0, "1")) +
            QConst.fQheadLS * g0(HeadLS.Q, 0, "\U00B5mol*m^-2*s^-1") +
            QConst.fQconsoleLS *
              g0(ConsoleLS.Q, 0, "\U00B5mol*m^-2*s^-1")
        }
      ),
      gfs3000 = list(
        desc = "GFS-3000 specific, for the top sensor.",
        deps = "",
        fn = function() { Meas.QambIn@"\U00B5mol*m^-2*s^-1" }
      ),
      gfs3000.gfs3000_light_bot = list(
        desc = "GFS-3000 specific, for the bottom sensor.",
        deps = "",
        fn = function() { Meas.QambInBot@"\U00B5mol*m^-2*s^-1" }
      )
    ),
    LeafQ.alpha = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { g0(LeafQ.alpha, LQConst.AbsAmbient) }
      ),
      ciras4 = list(
        desc = "",
        deps = "",
        fn = function() {
          # currently only light source weighted, not spectrum weighted
            ((g0(QConst.fQambIn) * g0(LQConst.AbsAmbient) +
               g0(QConst.fQambOut) * g0(LQConst.AbsAmbient) +
               g0(QConst.fQflr) * g0(LQConst.RedAbsFlr) +
               g0(QConst.fQheadLS) * g0(LQConst.RedAbsLED)) /
            (g0(QConst.fQambIn) + g0(QConst.fQambOut) +
               g0(QConst.fQflr) + g0(QConst.fQheadLS)))@.
        }
      ),
      li6400 = list(
        desc = "Li6400 specific, spectrum-weighted version.",
        deps = "",
        fn = function() {
          (g0(FlrLS.fblue, NA_real_, "1") * g0(LQConst.BlueAbsFlr, 0.85) +
            (1@"1" - g0(FlrLS.fblue, NA_real_, "1")) * g0(LQConst.RedAbsFlr, 0.85))@.
        }
      ),
      li6800 = list(
        desc = "Li6800 specific, spectrum-weighted version.",
        deps = "",
        fn = function() {
          un <- "\U00B5mol*m^-2*s^-1"
          rb <- pmax(g0(FlrLS.Qred, 0, un) + g0(FlrLS.Qmodavg, 0, un) +
                       g0(FlrLS.Qblue, 0, un),
                     set_units(0.1, un, mode = "standard"))
          flr <- ((g0(FlrLS.Qred, 0, un) + g0(FlrLS.Qmodavg, 0, un)) /
                    rb * LQConst.RedAbsFlr + g0(FlrLS.Qblue, 0, un) /
                    rb * LQConst.BlueAbsFlr)@.
          hls <- (g0(HeadLS.fred, 0, "1") * LQConst.RedAbsLED +
                    g0(HeadLS.fgreen, 0, "1") * LQConst.GreenAbsLED +
                    g0(HeadLS.fblue, 0, "1") * LQConst.BlueAbsLED +
                    g0(HeadLS.fwhite, 0, "1") * LQConst.WhiteAbsLED)@.
          ((QConst.fQambIn * LQConst.AbsAmbient +
              QConst.fQambOut * LQConst.AbsAmbient +
              QConst.fQflr * flr + QConst.fQheadLS * hls) /
              (QConst.fQambIn + QConst.fQambOut + QConst.fQflr +
                 QConst.fQheadLS))@.
        }
      )
    ),
    LeafQ.Conv = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { g0(LeafQ.Conv, g0(ConvAmbient, 0.19),
                             "J*\U00B5mol^-1")
        }
      ),
      ciras4 = list(
        desc = "",
        deps = "",
        fn = function() {
          # currently only light source weighted, not spectrum weighted
            ((g0(QConst.fQambIn) * g0(LQConst.ConvAmbient) +
               g0(QConst.fQambOut) * g0(LQConst.ConvAmbient) +
               g0(QConst.fQflr) * g0(LQConst.RedConvFlr) +
               g0(QConst.fQheadLS) * g0(LQConst.RedConvLED)) /
            (g0(QConst.fQambIn) + g0(QConst.fQambOut) +
               g0(QConst.fQflr) + g0(QConst.fQheadLS)))@"J*\U00B5mol^-1"
        }
      ),
      li6800 = list(
        desc = "Li6800 specific, accounts for the light spectrum.",
        deps = "",
        fn = function() {
          #NB: not all permutations tested!
          un <- "\U00B5mol*m^-2*s^-1"
          rb <- pmax(g0(FlrLS.Qred, 0, un) + g0(FlrLS.Qmodavg, 0, un) +
                       g0(FlrLS.Qblue, 0, un), 0.1@un)
          am <- QConst.fQambIn * LQConst.ConvAmbient + QConst.fQambOut *
            LQConst.ConvAmbient
          flr <- QConst.fQflr *
            ((g0(FlrLS.Qred, 0, un) + g0(FlrLS.Qmodavg, 0, un)) /
               rb * LQConst.RedConvFlr + g0(FlrLS.Qblue, 0, un) / rb *
               LQConst.BlueConvFlr)
          led <- QConst.fQheadLS *
            (g0(HeadLS.fred) * LQConst.RedConvLED +
               g0(HeadLS.fgreen) * LQConst.GreenConvLED +
               g0(HeadLS.fblue) * LQConst.BlueConvLED +
               g0(HeadLS.fwhite) * LQConst.WhiteConvLED)
          con <- QConst.fQconsoleLS *
            (g0(ConsoleLS.fred) * LQConst.RedConvLED +
               g0(ConsoleLS.fgreen) * LQConst.GreenConvLED +
               g0(ConsoleLS.fblue) * LQConst.BlueConvLED +
               g0(ConsoleLS.fwhite) * LQConst.WhiteConvLED)
          ((am + flr + led + con) /
              (QConst.fQambIn + QConst.fQambOut + QConst.fQflr +
                 QConst.fQheadLS + QConst.fQconsoleLS))@"J*\U00B5mol^-1"
        }
      )
    ),
    LeafQ.Qabs = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { LeafQ.Qin * LeafQ.alpha }
      )
    ),
    Dynamic.Adyn = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          if (exists("Dynamic.Crd") && exists("Dynamic.Csd") &&
              exists("Dynamic.alphaVc") && exists("Dynamic.dCsd_dt"))
            (Meas.Flow / Const.S * ((unity - Meas.H2Or) / unity) *
               (Dynamic.Crd - Dynamic.Csd - Meas.Pa * Dynamic.alphaVc *
                  Dynamic.dCsd_dt / (gasconstant * Meas.Tair@"K" *
                                       Meas.Flow)))@"\U00B5mol*m^-2*s^-1"
          else if(exists("Dynamic.Hr")) #gfs doesn't always use CO2buf?
            ((g0(Dynamic.Crd, Meas.CO2r) - Meas.CO2s *
               (unity - Dynamic.Hr) / (unity - Meas.H2Os)) * Meas.Flow /
               Const.S)@"\U00B5mol*m^-2*s^-1"
          else
            NA_real_@"\U00B5mol*m^-2*s^-1"
          }
      )
    ),
    Dynamic.Edyn = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          if (exists("Dynamic.Hs") && exists("Dynamic.Hr") &&
              exists("Dynamic.alphaVc") && exists("Dynamic.dHs_dt"))
            (Meas.Flow / Const.S * (unity / (unity - Dynamic.Hs)) *
               (Dynamic.Hs - Dynamic.Hr + Meas.Pa * Dynamic.alphaVh *
                  Dynamic.dHs_dt / (gasconstant * Meas.Tair@"K" *
                                      Meas.Flow)))@"mol*m^-2*s^-1"
          else if(exists("Dynamic.Hr"))
            (Meas.Flow * (Meas.H2Os - Dynamic.Hr) /
               (Const.S * (unity - Meas.H2Os)))@"mol*m^-2*s^-1"
          else
            NA_real_@"mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.gbw = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(GasEx.gbw, g0(Const.CustomBLC, 2, "mol*m^-2*s^-1"),
             "mol*m^-2*s^-1")
        }
      ),
      li6400 = list(
        desc = "Li6400 specific version. Does not account for fan speed.",
        deps = "",
        fn = function() {
          (Const.S@. * ChambConst.SslopeBLC +
             ChambConst.SoffsetBLC)@"mol*m^-2*s^-1"
        }
      ),
      li6800 = list(
        desc = "Li6800 specific version.",
        deps = "",
        fn = function() {
          (ifelse(substr(SysConst.Geometry, 1, 1) == "0",
                  { fpo <- (Meas.FanSpeed@"krpm" * Meas.Pa /
                    ChambConst.PoBLC)@.
                  mmS <- (pmax(pmin(Const.S, ChambConst.SmaxBLC),
                              ChambConst.SminBLC))@"cm^2"@.
                  ChambConst.CFaBLC + fpo * (ChambConst.CFbBLC +
                                              ChambConst.CFcBLC *
                                              mmS^2 + ChambConst.CFdBLC * mmS +
                                              ChambConst.CFeBLC * fpo)
                  }, ifelse(substr(SysConst.Geometry, 1, 1) == "1", 3,
                            as.numeric(Const.CustomBLC))))@"mol*m^-2*s^-1"
        }
      ),
      gfs3000 = list(
        desc = paste0("GFS-3000 specific version. Only valid for 8 cm^2 and ",
                      " the default chamber."),
        deps = "",
        fn = function() {
          # Based on a fit of data provided by the vendor for MEBA191A
          # Note the value is divided by 2 because we define gb as one-sided
          # Area correction is WIP and not yet here (there are issues)
          # FIXME: depend on chamber type?
          imp <- pmax(1, Meas.FanSpeed@steps@.)
          if (any(Const.S@"cm^2"@. != 8))
            warning("\n  No area correction for gbw is ",
                    "currently implemented.\n")
          ((-52.2029 * imp^2 + 1161.15 * imp - 538.845) / 2000)@"mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.VPcham = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { Meas.H2Os * (Meas.Pa + Meas.DeltaPcham) }
      )
    ),
    GasEx.SVPcham = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          # this eq is used by the 6400 and 6800
          0.61365 * exp(17.502 * GasEx.TairCnd@. /
                          (240.97 + GasEx.TairCnd@.))@"kPa"
        }
      ),
      GoffGratch1946 = list(
        desc = paste0("Based on the Goff-Gratch equation (used internally ",
                      "by the GFS-3000)."),
        deps = "",
        fn = function() {
          tr <- 373.16 / GasEx.TairCnd@"K"@.
          (10^(-7.90298 * (tr - 1) + 5.02808 * log10(tr) -
                 1.3816e-07 * (10^(11.344 * (1 - 1 / tr)) - 1) +
                 0.0081328 * (10^(-3.49149 * (tr - 1)) - 1) +
                 log10(1013.246) - 1))@"kPa"
        }
      ),
      Buck1981 = list(
        desc = paste0("This version uses the orginal (1981) Buck equation. ",
                      "It includes a correction for air pressure ",
                      "inside the leaf chamber."),
        deps = "",
        fn = function() {
          P <- (Meas.Pa + Meas.DeltaPcham)@"kPa"@"mbar"@.
          ((1.0007 + 3.46e-6 * P) * 0.61121 *
              exp(17.502 * GasEx.TairCnd@. / (240.97 + GasEx.TairCnd@.)))@"kPa"
        }
      ),
      Buck1996 = list(
        desc = paste0("This version uses the updated Buck (1996) equation. ",
                      "It includes a correction for air pressure ",
                      "inside the leaf chamber."),
        deps = "",
        fn = function() {
          tc <- GasEx.TairCnd@.
          P <- (Meas.Pa + Meas.DeltaPcham)@"kPa"@"mbar"@.
          EF <- 1 + 1e-4 * (7.2 + P * (0.0320 + 5.9e-6 * tc^2))
          (EF * (0.61121 * exp((18.678 - tc/234.5) * (tc/(257.14 + tc)))))@"kPa"
        }
      )
    ),
    GasEx.Rabs = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (LeafQ.Qin * LeafQ.Conv)@"W*m^-2" }
      ),
      li6400 = list(
        desc = paste0("Li6400 specific version that is corrected for",
                      "radiation leaving the chamber."),
        deps = "",
        fn = function() { ((Meas.QambIn * QConst.fQin +
                              Meas.QambOut * QConst.fQout) *
                             LeafQ.Conv)@"W*m^-2" }
      )
    ),
    GasEx.RHcham = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (GasEx.VPcham / GasEx.SVPcham)@"%" }
      )
    ),
    GasEx.Asty = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          (Meas.Flow * Leak.CorrFact *
             (Meas.CO2r - Meas.CO2s * (unity - Leak.CorrFact * Meas.H2Or) /
                (unity - Leak.CorrFact * Meas.H2Os)) /
             Const.S)@"\U00B5mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.Esty = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          (Meas.Flow * Leak.CorrFact * (Meas.H2Os - Meas.H2Or) /
             (Const.S * (unity - Leak.CorrFact * Meas.H2Os)))@"mmol*m^-2*s^-1"
        }
      )
    ),
    GasEx.Emm = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          ifelse(SysConst.UseDynamic@"mmol*m^-2*s^-1",
                 Dynamic.Edyn@"mmol*m^-2*s^-1",
                 GasEx.Esty@"mmol*m^-2*s^-1")
        }
      )
    ),
    GasEx.E = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { GasEx.Emm@"mol*m^-2*s^-1" }
      )
    ),
    GasEx.A = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          #default to steady state for GFS:
          ifelse(SysConst.UseDynamic &
                   SysObs.Instrument != "GFS3000",
                 Dynamic.Adyn@"\U00B5mol*m^-2*s^-1",
                 GasEx.Asty@"\U00B5mol*m^-2*s^-1")@"\U00B5mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.Ca = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { Meas.CO2s@"\U00B5mol*mol^-1" }
      ),
      li6800 = list(
        desc = "Li6800 specific version that includes a leak correction.",
        deps = "",
        fn = function() {
          # note 6800 firmware < 2.1.11 had a small bug
          (Meas.CO2s - ifelse((Leak.CorrFact > 1)@"\U00B5mol*mol^-1",
                              (GasEx.A * Const.S / Leak.Fan)@"\U00B5mol*mol^-1",
                              0))@"\U00B5mol*mol^-1"
        }
      )
    ),
    GasEx.pCa = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (GasEx.Ca * (Meas.Pa + Meas.DeltaPcham))@"Pa" }
      )
    ),
    GasEx.LatHFlux = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          # might define 44100 as a constant in future
          44100@"W*s*mol^-1" * -GasEx.E
        }
      )
    ),
    GasEx.SenHFlux = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          #29.3 and 0.92 may be variable, and could be settable constants?
          (2 * 29.3@"J*mol^-1*K^-1" * GasEx.gbw * 0.92 *
             (GasEx.TairCnd@"K" - GasEx.TleafCnd@"K"))@"W*m^-2"
        }
      )
    ),
    GasEx.NetTherm = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          (2 * 0.95 * 5.67e-08@"W*K^-4*m^-2" *
             ((Meas.Tair@"K" + LTConst.DeltaTw)^4 -
                (GasEx.TleafCnd@"K")^4))@"W*m^-2"
        }
      )
    ),
    GasEx.EBSum = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          GasEx.Rabs + GasEx.NetTherm + GasEx.LatHFlux + GasEx.SenHFlux
        }
      )
    ),
    GasEx.TleafEB = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          epssigma <-  0.95 * 5.67e-08@"W*K^-4*m^-2"
          Ta <- GasEx.TairCnd@"K"
          # Tair (not TairCnd) for Tw!
          Tw <- Meas.Tair@"K" + LTConst.DeltaTw
          (Ta + (GasEx.Rabs + 2 * epssigma * (Tw^4 - Ta^4) + GasEx.LatHFlux) /
              (2 * 29.3@"J*mol^-1*K^-1" * 0.92 * GasEx.gbw +
                 8 * epssigma * Ta^3))@"degC"
        }
      )
    ),
    GasEx.TairCnd = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { Meas.Tair }
      ),
      li6400 = list(
        desc = paste0("The air temperature for the Li6400 is estimated ",
                      "as the average between the leaf and air temperature ",
                      "sensors or measured with the leaf thermocouple ",
                      "if it is not touching the leaf in energy-balance mode."),
        deps = "",
        fn = function() { (ifelse(LTConst.fTEB@., Meas.Tleaf,
                                  (Meas.Tleaf + Meas.Tair) / 2))@"degC" }
      )
    ),
    GasEx.TleafCnd = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          LTConst.fT1 * g0(Meas.Tleaf, -999.99, "degC", TRUE) +
            LTConst.fT2 * g0(Meas.Tleaf2, -999.99, "degC", TRUE) +
            LTConst.fTEB * g0(GasEx.TleafEB, -999.99, "degC", TRUE)
        }
      )
    ),
    GasEx.SVPleaf = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          0.61365 * exp(17.502 * GasEx.TleafCnd@. /
                          (240.97 + GasEx.TleafCnd@.))@"kPa"
        }
      ),
      GoffGratch1946 = list(
        desc = paste0("Based on the Goff-Gratch equation (used internally ",
                      "by the GFS-3000)."),
        deps = "",
        fn = function() {
          tr <- 373.16 / GasEx.TleafCnd@"K"@.
          (10^(-7.90298 * (tr - 1) + 5.02808 * log10(tr) -
                 1.3816e-07 * (10^(11.344 * (1 - 1 / tr)) - 1) +
                 0.0081328 * (10^(-3.49149 * (tr - 1)) - 1) +
                 log10(1013.246) - 1))@"kPa"
        }
      ),
      Buck1981 = list(
        desc = paste0("This version uses the orginal (1981) Buck equation. ",
                      "It includes a correction for air pressure ",
                      "inside the leaf chamber."),
        deps = "",
        fn = function() {
          ((1.0007 + 3.46e-6 * (Meas.Pa + Meas.DeltaPcham)@"kPa"@"mbar"@.) *
             0.61121 * exp(17.502 * GasEx.TleafCnd@. /
                             (240.97 + GasEx.TleafCnd@.)))@"kPa"
        }
      ),
      Buck1996 = list(
        desc = paste0("This version uses the updated Buck (1996) equation. ",
                      "It includes a correction for air pressure ",
                      "inside the leaf chamber."),
        deps = "",
        fn = function() {
          tc <- GasEx.TleafCnd@.
          P <- (Meas.Pa + Meas.DeltaPcham)@"kPa"@"mbar"@.
          EF <- 1 + 1e-4 * (7.2 + P * (0.0320 + 5.9e-6 * tc^2))
          (EF * (0.61121 * exp((18.678 - tc / 234.5) *
                                 (tc / (257.14 + tc)))))@"kPa"
        }
      )
    ),
    GasEx.VPDleaf = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          (GasEx.SVPleaf * g0(Const.RHi, 1, "1") -
             Meas.H2Os * (Meas.Pa + Meas.DeltaPcham))@"kPa"
        }
      )

    ),
    GasEx.gtw = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          wi <- GasEx.SVPleaf * g0(Const.RHi, 1, "1") /
            (Meas.Pa + Meas.DeltaPcham)
          (GasEx.E  * (unity - (wi + Meas.H2Os) / 2) /
              (wi - Meas.H2Os))@"mol*m^-2*s^-1"
        }
      ),
      cuticular_conductance = list(
        desc = paste0("This version is based on Marquez et al. 2021, and ",
                      "accounts for cuticular conductance. It has been ",
                      "modified to account for differences in adaxial ",
                      "and abaxial conductance."),
        deps = "",
        fn = function() {
          # Marquez, Stuart-Williams and Farquhar 2021 Eq. 9 and 10
          # but split into an adaxial and abaxial part
          # also note licor gbw is ONE-SIDED and is doubled here.
          wi <- GasEx.SVPleaf * g0(Const.RHi, 1, "1") /
                   (Meas.Pa + Meas.DeltaPcham)
          ws <- (GasEx.E / (2 * GasEx.gbw) * (unity - Meas.H2Os / 2) +
                   Meas.H2Os) / (unity + GasEx.E / (4 * GasEx.gbw))
          Es <- GasEx.E - Const.gcw * (wi - ws)
          gtw_up <- (1 / ((wi - ws) * (1 + Const.K) /
                            (GasEx.E - Es * (wi + ws) / 2) +
                            1 / (GasEx.gbw)))
          gtw_low <- (Const.K / ((wi - ws) * (1 + Const.K) /
                                   (GasEx.E - Es * (wi + ws) / 2) +
                                   Const.K / (GasEx.gbw)))
          (gtw_up + gtw_low)@"mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.gsw = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          rt <- 1 / GasEx.gtw
          rb <- 1 / GasEx.gbw
          (2 / (rt - rb + sign(rt) * sqrt((rt - rb) ^ 2 + 4 * Const.K /
                                            (Const.K + 1) ^ 2 * (2 * rt * rb -
                                                                   rb ^ 2)))
          )@"mol*m^-2*s^-1"
        }
      ),
      cuticular_conductance = list(
        desc = paste0("This version is based on Marquez et al. 2021, and ",
                      "accounts for cuticular conductance. It has been ",
                      "modified to account for differences in adaxial ",
                      "and abaxial conductance."),
        deps = "",
        fn = function() {
          rt <- 1 / GasEx.gtw
          rb <- 1 / GasEx.gbw
          k <- pmax(Const.K, 1e-12) # avoid /0
          Q <-  2 * k * rb * (rb - 2 * rt) / (1 + k)^2
          ((rt - rb - Q * Const.gcw - sign(rt) * sqrt(-2 * Q + (rb - rt)^2)) /
              Q)@"mol*m^-2*s^-1"
        }
      )
    ),
    GasEx.gtc = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          gs <- GasEx.gsw / 1.6
          gb <- GasEx.gbw / (1.6)^(2 / 3)
          gs * gb * (1 / ((Const.K + 1) * gb + gs ) + Const.K /
                       ((Const.K + 1) * gb  + Const.K * gs))
        }
      ),
      cuticular_conductance = list(
        desc = paste0("This version is based on ideas in Marquez et al. 2021, ",
                      "and accounts for cuticular conductance and differences ",
                      "in adaxial and abaxial conductance."),
        deps = "",
        fn = function() {
          # add a cuticular conductance in the spirit
          # of Marquez et al. 2021
          gl <- GasEx.gsw / 1.6 + Const.gcc
          gb <- GasEx.gbw / (1.6)^(2 / 3)
          k <- pmax(Const.K, 1e-12)
          1 / (1 / gb + (1 + k) / gl) + k / (k / gb + (1 + k) / gl)
        }
      )
    ),
    GasEx.Cs = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          # Based on Cs-Mix of Marquez thesis and Note 6 Marquez et al. 2023
         ((GasEx.Ca * (2 * GasEx.gbw / (1.6^(2 / 3)) - GasEx.E / 2) -
                   GasEx.A) / (2 * GasEx.gbw / (1.6^(2 / 3)) +
                                 GasEx.E / 2))@"\U00B5mol*mol^-1"
        }
      )
    ),
    GasEx.Ci = list(
      default = list(
        desc = "",
        deps = "",
        fn =  function() {
          (((GasEx.gtc - GasEx.E / 2) * GasEx.Ca - GasEx.A) /
             (GasEx.gtc + GasEx.E / 2))@"\U00B5mol*mol^-1"
        }
      ),
      cuticular_conductance = list(
        desc = paste0("This version is based on Marquez 2021, and ",
                      "Marquez et al. 2023. The resulting Ci is somewhere ",
                      "between adaxial and abaxial values."),
        deps = "",
        fn = function() {
          glc <- Const.gcc + GasEx.gsw / 1.6
          # note licor gbw is ONE-SIDED and is doubled here.
          wi <- GasEx.SVPleaf * g0(Const.RHi, 1, "1") /
            (Meas.Pa + Meas.DeltaPcham)
          ws <- (GasEx.E / (2 * GasEx.gbw) * (unity - Meas.H2Os / 2) +
                   Meas.H2Os) / (unity + GasEx.E / (4 * GasEx.gbw))
          Es <- GasEx.E - Const.gcw * (wi - ws)
          ((GasEx.Cs * (glc - Es / 2) - GasEx.A) /
              (glc + Es / 2))@"\U00B5mol*mol^-1"
        }
      )
    ),
    GasEx.Ci_Ca = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {GasEx.Ci / GasEx.Ca}
      )
    ),
    GasEx.pCi = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (GasEx.Ci * (Meas.Pa + Meas.DeltaPcham))@"Pa" }
      )
    ),
    Meas.H2Oa = list(
      raw.li6800 = list(
        desc = paste0("Li6800 specific version that uses raw absorption",
                      "values and accounts for the effect of oxygen, ",
                      "temperature and pressure."),
        deps = "match",
        fn = function() {
          #FIXME: error on missing cal? or Raw?
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          ab <- span_abs(hc, "h2o", "a", Meas.Pa, Raw.H2OaAbs)
          (abs2frac(hc, "h2o", "a", ab) * Status.Ts@"K" *
             psiH2O(hc, Const.Oxygen@"%"))@"mmol*mol^-1"
        }
      ),
      O2_correction.li6800 = list(
        desc = paste0("Li6800 specific version that back-calculates mol ",
                      "fractions and takes the effect of a possible change in ",
                      "oxygen concentration into account."),
        deps = "match",
        fn = function() {
          #FIXME: warning on missing cal?
          (Raw.H2Oa * psiH2O(get_cal(SysConst.UserCal, SysConst.FactCal),
                              Const.Oxygen@"%"))@"mmol*mol^-1"
        }
      ),
      O2_correction.gfs3000 = list(
        desc = paste0("GFS-3000 specific version that back-calculates mol ",
                      "fractions and takes the effect of a possible change in ",
                      "oxygen concentration into account."),
        deps = "match",
        fn = function() {
          # base factor on H2Oa instead? Effect is minor...
          (Raw.H2Oa *
             gfs_o2_factor(Raw.H2Or, Const.Oxygen@"%"))@"mmol*mol^-1"
        }
      )
    ),
    Meas.H2Or = list(
      raw.li6400 = list(
        desc = "",
        deps = "",
        fn = function() {
          warning("Calculating mole fractions from raw values not yet ",
                  "implemented.")
          Meas.H2Or
          }
      ),
      raw.ciras4 = list(
        desc = "",
        deps = "",
        fn = function() {
          warning("Calculating mole fractions from raw values not possible.")
          Meas.H2Or
        }
      ),
      raw.gfs3000 = list(
        desc = "",
        deps = "match",
        fn = function() {
          warning("Calculating mole fractions from raw values not possible.")
          Meas.H2Or
        }
      ),
      raw.li6800 = list(
        desc = paste0("Li6800 specific version that uses raw absorption",
                      "values and accounts for the effect of oxygen, ",
                      "temperature and pressure."),
        deps = "match",
        fn = function() {
          #FIXME: error on missing cal? or Raw?
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          ab <- span_abs(hc, "h2o", "b", Meas.Pa, Raw.H2OrAbs)
          (abs2frac(hc, "h2o", "b", ab) * Status.Tr@"K" *
              psiH2O(hc, Const.Oxygen@"%"))@"mmol*mol^-1"
        }
      ),
      O2_correction.li6400 = list(
        desc = "",
        deps = "",
        fn = function() {
          warning("O2 corrections for Li6400 are not yet implemented.")
          Meas.H2Or
        }
      ),
      O2_correction.ciras4 = list(
        desc = "",
        deps = "",
        fn = function() {
          warning("O2 corrections for CIRAS4 are not yet implemented.")
          Meas.H2Or
        }
      ),
      O2_correction.li6800 = list(
        desc = paste0("Li6800 specific version that back-calculates mol ",
                      "fractions and takes the effect of a possible change in ",
                      "oxygen concentration into account."),
        deps = "match",
        fn = function() {
          #FIXME: warning on missing cal
          (Raw.H2Or * psiH2O(get_cal(SysConst.UserCal, SysConst.FactCal),
                              Const.Oxygen@"%"))@"mmol*mol^-1"
        }
      ),
      O2_correction.gfs3000 = list(
        desc = paste0("GFS-3000 specific version that back-calculates mol ",
                      "fractions and takes the effect of a possible change in ",
                      "oxygen concentration into account."),
        deps = "match",
        fn = function() {
          (Raw.H2Or * gfs_o2_factor(Raw.H2Or,
                                     Const.Oxygen@"%"))@"mmol*mol^-1"

        }
      )
    ),
    Meas.CO2a = list(
      raw.li6800 = list(
        desc = paste0("Calculate sample CO2 mol fraction from raw ",
                      "absorption values, acCounting for the effect of ",
                      "oxygen, temperature and pressure."),
        deps = "match",
        fn = function() {
          #FIXME: error on missing cal? or Raw?
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          p <- psiCO2(hc, Const.Oxygen, Meas.H2Oa)
          ab <- span_abs(hc, "co2", "a", Meas.Pa, Raw.CO2aAbs)
          (abs2frac(hc, "co2", "a", ab / p) *
              Status.Ts@"K" * p)@"\U00B5mol*mol^-1"

        }
      ),
      O2_correction.li6800 = list(
        desc = paste0("Recalculate sample CO2 mol fraction for the 6800, ",
                      "taking the effect of a possible change in oxygen ",
                      "concentration into acCount."),
        deps = "match",
        fn = function() {
          #FIXME: warning on missing cal
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          p <- psiCO2(hc, Const.Oxygen, Meas.H2Oa)
          (abs2frac(hc, "co2", "a", Raw.CO2aAbsP / p) * Status.Ts@"K" *
              p)@"\U00B5mol*mol^-1"
        }
      )
    ),
    Meas.CO2r = list(
      raw.li6800 = list(
        desc = paste0("Li6800 specific version that uses raw absorption",
                       "values and accounts for the effect of oxygen, ",
                       "temperature and pressure."),
        deps = "match",
        fn = function() {
          #FIXME: error on missing cal? or Raw?
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          p <- psiCO2(hc, Const.Oxygen, Meas.H2Or)
          ab <- span_abs(hc, "co2", "b", Meas.Pa, Raw.CO2rAbs)
          (abs2frac(hc, "co2", "b", ab / p) *
              Status.Tr@"K" * p)@"\U00B5mol*mol^-1"

        }
      ),
      O2_correction.li6800 = list(
        desc = paste0("Li6800 specific version that uses raw absorption",
                      "values and accounts for the effect of oxygen, ",
                      "temperature and pressure."),
        deps = "match",
        fn = function() {
          #FIXME: warning on missing cal
          hc <- get_cal(SysConst.UserCal, SysConst.FactCal)
          p <- psiCO2(hc, Const.Oxygen, Meas.H2Or)
          (abs2frac(hc, "co2", "b", Raw.CO2rAbsP / p) * Status.Tr@"K" *
              p)@"\U00B5mol*mol^-1"
        }
      )
    ),
    Meas.CO2s = list(
      match.li6800 = list(
        desc = paste0("Li6800 specific version that recalculates sample CO2 ",
                      "mol fraction using the match correction."),
        deps = "",
        fn = function() {
          #despite what the manual suggests, the match function is only used
          #for range matching. But there is no flag to indicate range is active.
          #for now we rely on the stored offsets:
          Meas.CO2a + MchStatus.MatchCO2
        }
      ),
      match.gfs3000 = list(
        desc = paste0("GFS-3000 specific version that recalculates sample ",
                      "CO2 mol fraction using the GFS-3000 match correction."),
        deps = "",
        fn =  function() {
          (Meas.CO2a@"\U00B5mol*mol^-1" +
             MchEvent.CO2match@"\U00B5mol*mol^-1")
        }
      )
    ),
    Meas.H2Os = list(
      match.li6800 = list(
        desc = paste0("Li6800 specific version that recalculates sample H2O ",
                      "mol fraction using a match correction function."),
        deps = "",
        fn = function() {
          #despite what the manual suggests, the match function is only used
          #for range matching. But there is no flag to indicate range is active.
          #for now we rely on the stored offsets:
          Meas.H2Oa + MchStatus.MatchH2O
        }
      ),
      match.gfs3000 = list(
        desc =  paste0("GFS-3000 specific version that recalculates sample ",
                       "H2O mol fraction using the GFS-3000 match correction."),
        deps = "",
        fn = function() {
          (Meas.H2Oa@"mmol*mol^-1" + MchEvent.H2Omatch@"mmol*mol^-1")
        }
      )
    ),
    FLR.Fv_Fm = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          Fv_Fm <- 1 - g0(FLR.Fo, NA_real_) / g0(FLR.Fm, NA_real_)
        }
      )
    ),
    FLR.Fvp_Fmp = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { 1 - g0(FLR.Fop, NA_real_)/g0(FLR.Fmp, NA_real_) }
      )
    ),
    FLR.phiPS2 = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { 1 - g0(FLR.Fs, NA_real_)/g0(FLR.Fmp, NA_real_) }
      )
    ),
    FLR.phiCO2 = list(
      default = list(
        desc = "Quantum efficency of the carboxylation rate.",
        deps = "",
        fn = function() {
          ((GasEx.A - g0(FLR.Adark, NA_real_, "\U00B5mol*m^-2*s^-1")) /
            FLR.QabsFs)@"\U00B5mol*mol^-1" }
      )
    ),
    FLR.QabsFs = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(FLR.QinFs, LeafQ.Qin, "\U00B5mol*m^-2*s^-1") * LeafQ.alpha
        }
      )
    ),
    FLR.ETR = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(FLR.phiPS2, NA_real_) * g0(Const.fPS2, 0.5) * FLR.QabsFs
          }
      )
    ),
    FLR.NPQ = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (g0(FLR.Fm, NA_real_) - g0(FLR.Fmp, NA_real_)) /
            g0(FLR.Fmp, NA_real_)
          }
      )
    ),
    FLR.qN = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (g0(FLR.Fm, NA_real_) - g0(FLR.Fmp, NA_real_)) /
            (g0(FLR.Fm, NA_real_) - g0(FLR.Fop, NA_real_))
        }
      )
    ),
    FLR.qNFo = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (g0(FLR.Fm, NA_real_) - g0(FLR.Fmp, NA_real_)) /
            (g0(FLR.Fm, NA_real_) - g0(FLR.Fo, NA_real_))
        }
      )
    ),
    FLR.qP = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (g0(FLR.Fmp, NA_real_) - g0(FLR.Fs, NA_real_)) /
            (g0(FLR.Fmp, NA_real_) - g0(FLR.Fop, NA_real_))
        }
      )
    ),
    FLR.qPFo = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { (g0(FLR.Fmp, NA_real_) - g0(FLR.Fs, NA_real_)) /
            (g0(FLR.Fmp, NA_real_) - g0(FLR.Fo, NA_real_))
        }
      )
    ),
    FLR.qL = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { g0(FLR.qP, NA_real_) * g0(FLR.Fop, NA_real_) /
            g0(FLR.Fs, NA_real_)
        }
      )
    ),
    FLR.Fop = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          fopalt <- g0(FLR.FopAlt, NA_real_)
          fmin <- g0(FLR.Fmin, na_replace = TRUE) * rep(1, length(fopalt))
          fmin[fmin == 0] <- g0(FLR.FopAlt, NA_real_)[fmin == 0]
        }
      )
    ),
    FLR.FopAlt = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {
          g0(FLR.Fo, NA_real_) /
            (g0(FLR.Fv_Fm, NA_real_) + g0(FLR.Fo, NA_real_) /
               g0(FLR.Fmp, NA_real_))
        }
      )
    ),
    FLR.PhiQin_4 = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() {FLR.phiPS2 * LeafQ.Qin / 4}
      )
    ),
    FLR.Cc = list(
      gm_fluorescence = list(
        desc = paste0("Chloroplast CO2 mol fractions calculated from electron ",
                      "transport rates estimated by chlorophyll fluorescence.",
                      "This estimate requires calibration of FLR.ETR and ",
                      "accurate estimates for Const.RL and Const.GammaStar."),
        deps = "",
        fn = function() {
          GP <- GasEx.A + g0(Const.RL, NA_real_, "\U00B5mol*m^-2*s^-1")
          JF <- g0(FLR.ETR, NA_real_, "\U00B5mol*m^-2*s^-1")
          g0(Const.GammaStar, NA_real_) * (JF + 8 * GP) / (JF - 4 * GP)
        }
      )
    ),
    FLR.gm = list(
      gm_fluorescence = list(
        desc = paste0("Mesophyll conductance calculated from electron ",
                      "transport rates estimated by chlorophyll fluorescence.",
                      "This estimate requires calibration of FLR.ETR and ",
                      "accurate estimates for Const.RL and Const.GammaStar."),
        deps = "",
        fn = function() {
          GasEx.A / (GasEx.Ci - FLR.Cc)
        }
      )
    ),
    d13C.A_pCa = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() { GasEx.A / GasEx.pCa }
      )
    ),
    d13C.ap = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          (d13CConst.ab * (GasEx.Ca - GasEx.Cs) + d13CConst.as *
             (GasEx.Cs - GasEx.Ci)) / (GasEx.Ca - GasEx.Ci) }
      )
    ),
    d13C.t = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          ((1@"1" + d13C.ap) * GasEx.E)/(2 * GasEx.gtc)
        }
      )
    ),
    d13C.ep = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          d13CConst.e + g0(d13CMeas.delta13CO2r, 0, "permille") -
            g0(d13CConst.delta13CO2g, 0, "permille")
          }
      ),
      d13C_e_Busch2020.d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          d13CConst.e + (d13CMeas.delta13CO2s - d13C.Deltao) -
            (d13CConst.delta13CO2g - d13CConst.Deltag)
        }
      )
    ),
    d13C.xi = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          (Meas.CO2r / (1@"1" - Meas.H2Or)) /
                (Meas.CO2r / (1@"1"- Meas.H2Or) -
                   Meas.CO2s / (1@"1" - Meas.H2Os))
        }
      )
    ),
    d13C.Deltao = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          d13CO2s <- g0(d13CMeas.delta13CO2s, NA_real_, "permille")
          d13CO2r <- g0(d13CMeas.delta13CO2r, NA_real_, "permille")
          d13C.xi * (d13CO2s - d13CO2r) /
            (1@"1" + d13CO2s - d13C.xi * (d13CO2s - d13CO2r))

        }
      )
    ),
    d13C.Deltai = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          GS <- g0(Const.GammaStar, NA_real_, "\U00B5mol*mol^-1")
          RL <- g0(Const.RL, NA_real_, "\U00B5mol*m^-2*s^-1")
          t <- d13C.t@.
          one <- 1@"1"
          #Farquhar and Cernusak 2012, Eq 11
          1 / ( 1 - t) * d13C.ap * (GasEx.Ca - GasEx.Ci) / GasEx.Ca +
            (1 + t)/(1 - t) * (d13CConst.b * GasEx.Ci / GasEx.Ca -
                                 (one + d13CConst.b) / (one + d13C.ep) *
                                 d13C.ep * RL / (GasEx.A + RL) *
                                 (GasEx.Ci - GS) / GasEx.Ca -
                                 (one + d13CConst.b) / (one + d13CConst.f) *
                                 d13CConst.f * (GS / GasEx.Ca))

        }
      ),
      d13C_dis.d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          GS <- g0(Const.GammaStar, NA_real_, "\U00B5mol*mol^-1")
          RL <- g0(Const.RL, NA_real_, "\U00B5mol*m^-2*s^-1")
          t <- d13C.t@.
          one <- 1@"1"
          wh <- 0@"permille"
          #Busch et al 2020, S24, with added wh
          aR <- one + RL / GasEx.A * d13C.ep / (one + d13C.ep)
          1 / ( 1 - t) * d13C.ap * (GasEx.Ca - GasEx.Ci) / GasEx.Ca +
            (1 + t)/(1 - t) *
            (d13CConst.b * GasEx.Ci / GasEx.Ca -
               (one + d13CConst.b) / ((one + d13C.ep) * aR) *
               d13C.ep * RL / GasEx.A * GasEx.Ci / GasEx.Ca -
               (one + d13CConst.b) / ((one + d13CConst.f) * aR) *
               (d13CConst.f - wh) * (GS / GasEx.Ca))

        }
      )
    ),
    d13C.DeltaiDeltao = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() { (d13C.Deltai - d13C.Deltao) }
      )
    ),
    d13C.gm = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          RL <- g0(Const.RL, NA_real_, "\U00B5mol*m^-2*s^-1")
          t <- d13C.t@.
          one <- 1@"1"
          rm <- ((1 - t) / (1 + t) * d13C.DeltaiDeltao * GasEx.pCa /
                   (GasEx.A * (d13CConst.b - d13CConst.am -
                                 (one + d13CConst.b) /
                                 (one + d13C.ep) * d13C.ep *
                                 (RL / (GasEx.A + RL)))))
          (1 / rm)@"mol*m^-2*s^-1*bar^-1"
        }
      ),
      d13C_dis.d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          RL <- g0(Const.RL, NA_real_, "\U00B5mol*m^-2*s^-1")
          t <- d13C.t@.
          one <- 1@"1"
          #Busch et al 2020, S26
          aR <- one + RL / GasEx.A * d13C.ep / (one + d13C.ep)
          rm <- ((1 - t) / (1 + t) * d13C.DeltaiDeltao * GasEx.pCa /
                   (GasEx.A * (d13CConst.b - d13CConst.am -
                                 (one + d13CConst.b) /
                                 ((one + d13C.ep) * aR) * d13C.ep *
                                 (RL / GasEx.A))))
          (1 / rm)@"mol*m^-2*s^-1*bar^-1"
        }
      )
    ),
    d13C.Cc = list(
      d13C = list(
        desc = "",
        deps = "",
        fn = function() {
          GasEx.Ci - GasEx.A / (d13C.gm * Meas.Pa)
        }
      )
    ),
    gasanalyzer.UseEqUnits = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { TRUE }
      )
    ),
    SysObs.Instrument = list(
      default = list(
        desc = "",
        deps = "",
        fn = function() { g0(SysObs.Instrument, NA_character_) }
      ),
      li6400 = list(
        desc = "",
        deps = "",
        fn = function() {
          inst <- g0(SysObs.Instrument, NA_character_)
          if (!isTRUE(all(SysObs.Instrument == "Li6400")))
            warning("\n  Applying Li6400 specific calculations to ",
                    "rows not measured by an Li6400.\n")
          inst }
      ),
      li6800 = list(
        desc = "",
        deps = "",
        fn = function() {
          inst <- g0(SysObs.Instrument, NA_character_)
          if (!isTRUE(all(SysObs.Instrument == "Li6800")))
            warning("\n  Applying Li6800 specific calculations to ",
                    "rows not measured by an Li6800.\n")
          inst }
      ),
      gfs3000 = list(
        desc = "",
        deps = "",
        fn = function() {
          inst <- g0(SysObs.Instrument, NA_character_)
          if (!isTRUE(all(SysObs.Instrument == "GFS3000")))
            warning("\n  Applying GFS3000 specific calculations to ",
                    "rows not measured by a GFS3000.\n")
          inst }
      ),
      ciras4 = list(
        desc = "",
        deps = "",
        fn = function() {
          inst <- g0(SysObs.Instrument, NA_character_)
          if (!isTRUE(all(SysObs.Instrument == "CIRAS4")))
            warning("\n  Applying CIRAS4 specific calculations to ",
                    "rows not measured by a CIRAS4.\n")
          inst }
      )
    )
  ))
  out
}






