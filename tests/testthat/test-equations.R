test_that("check creating equations", {

  tb <- tibble(
    Const.S = 1,
    GasEx.A = 1, GasEx.Ci = 500, GasEx.gsw = 1,
    Const.RL = 1, Const.GammaStar = 1,
    FLR.ETR = 1,
    Meas.Pa = 100, Meas.DeltaPcham = 0, Meas.Tair = 25, Meas.Tleaf = 25,
    Meas.QambIn = 100, Meas.QambInBot = 100, Meas.QambOut = 10,
    ChambConst.SslopeBLC = 1, ChambConst.SoffsetBLC = 0,
    Meas.FanSpeed = 10,
    GasEx.TleafCnd = 25, GasEx.TairCnd = 25,
    QConst.fQambIn = 1, LQConst.AbsAmbient = 0.5, LQConst.ConvAmbient = 0.5,
    QConst.fQambOut = 1, QConst.fQin = 1, QConst.fQout = 1,
    QConst.fQflr = 1, LQConst.RedAbsFlr = 0.5, LQConst.RedConvFlr = 0.5,
    LQConst.BlueAbsFlr = 0.5, FlrLS.fblue = 0.1,
    QConst.fQheadLS = 1, LQConst.RedAbsLED = 0.5, LQConst.RedConvLED = 0.5)

  expect_error(create_equations(""))
  # not sure how to test for 2 warnings, so I combine help with default:
  expect_warning(create_equations(c("default", "help")))

  nms <- names(create_equations("gm_fluorescence"))
  expect_setequal(nms , c("FLR.Cc", "FLR.gm", "gasanalyzer.UseFlags"))
  nms <- names(create_equations("Buck1981"))
  expect_setequal(nms , c("GasEx.SVPleaf", "GasEx.SVPcham",
                          "gasanalyzer.UseFlags"))
  nms <- names(create_equations("Buck1996"))
  expect_setequal(nms , c("GasEx.SVPleaf", "GasEx.SVPcham",
                          "gasanalyzer.UseFlags"))
  nms <- names(create_equations("GoffGratch1946"))
  expect_setequal(nms , c("GasEx.SVPleaf", "GasEx.SVPcham",
                          "gasanalyzer.UseFlags"))
  nms <- names(create_equations("ciras4"))
  expect_setequal(nms , c("SysObs.Instrument", "LeafQ.alpha",
                          "LeafQ.Conv", "gasanalyzer.UseFlags"))
  nms <- names(create_equations(c("gfs3000","gfs3000_light_bot")))
  expect_setequal(nms , c("SysObs.Instrument", "GasEx.gbw", "LeafQ.Qin",
                          "gasanalyzer.UseFlags"))
  nms <- names(create_equations("li6400"))
  expect_setequal(nms , c("SysObs.Instrument", "GasEx.gbw", "GasEx.Rabs",
                          "GasEx.TairCnd", "LeafQ.alpha",
                          "gasanalyzer.UseFlags"))
  nms <- names(create_equations("li6400"))
  expect_setequal(nms , c("SysObs.Instrument", "GasEx.gbw", "GasEx.Rabs",
                  "GasEx.TairCnd", "LeafQ.alpha", "gasanalyzer.UseFlags"))
  nms1 <- names(create_equations(c("li6800" , "raw")))
  nms2 <- names(create_equations(c("li6800" , "O2_correction")))
  expected <-  c("SysObs.Instrument", "GasEx.gbw", "LeafQ.Qin", "LeafQ.alpha",
                 "LeafQ.Conv", "Leak.Fan", "Meas.H2Oa", "Meas.H2Or",
                 "Leak.CorrFact", "Meas.CO2a", "Meas.CO2r", "Meas.H2Os",
                 "Meas.CO2s", "GasEx.Ca", "gasanalyzer.UseFlags")
  expect_setequal(nms1, expected)
  expect_setequal(nms2, expected)
  expect_error(create_equations(c("li6800", "O2_correction","raw")))
  nms <- names(create_equations(c("d13C", "d13C_dis", "d13C_e_Busch2020")))
  expect_setequal(nms , c("d13C.xi", "d13C.ap", "d13C.ep", "d13C.t",
                         "d13C.Deltai", "d13C.Deltao", "d13C.DeltaiDeltao",
                         "d13C.A_pCa", "d13C.gm", "d13C.Cc",
                         "gasanalyzer.UseFlags"))
  Eqs <- create_equations("default")

  examplefile <- system.file("extdata", "lowo2", package = "gasanalyzer")
  xlsxfile <- paste0(examplefile, ".xlsx")
  #FIXME: I should disable the warning
  xlsxdata <- suppressWarnings(read_6800_xlsx(xlsxfile))
  xlsxdata_norecalc <- read_6800_xlsx(xlsxfile, recalculate = FALSE)
  xlsx_eqs <- read_6800_equations(xlsxfile)
  txtdata <- read_6800_txt(examplefile)
  # xlsx stores zeroes for calced values, so only check constants and sys vals
  sys_cols <- names(xlsxdata)[grep("^Sys|Const" , names(xlsxdata))]
  # 1) For some reason, Li6800 Leak.Fan data differ slightly from the xlsx
  # 2) NA or NaN or Inf FLR variables get 0 in the txt
  check_cols <- names(xlsxdata)[!names(xlsxdata) %in%
                                  c("gasanalyzer.Equations",
                                    "Leak.Fan", "FLR.Fv_Fm", "FLR.FopAlt",
                                    "FLR.Fop", "FLR.Fvp_Fmp", "FLR.qP",
                                    "FLR.qN", "FLR.qNFo", "FLR.qL",
                                    "FLR.1_qL")]

  expect_equal(xlsx_eqs, xlsxdata$gasanalyzer.Equations[[1]],
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(xlsxdata[sys_cols], xlsxdata_norecalc[sys_cols],
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(xlsxdata[check_cols], txtdata[check_cols],
               tolerance = 1e-6, ignore_attr = TRUE)
})
