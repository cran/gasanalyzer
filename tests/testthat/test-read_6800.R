test_that("check read_6800_* methods", {
  options("gasanalyzer.calibration.warning" = FALSE)
  examplefile <- system.file("extdata//lowo2", package = "gasanalyzer")
  xlsxfile <- paste0(examplefile, ".xlsx")
  xlsxdata <- read_6800_xlsx(xlsxfile)
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

  expect_equal(xlsx_eqs, xlsxdata$gasanalyzer.Equations[[1]] )
  expect_equal(xlsxdata[sys_cols], xlsxdata_norecalc[sys_cols],
               tolerance = 1e-6)
  expect_equal(xlsxdata[check_cols], txtdata[check_cols],
               tolerance = 1e-6)
})
