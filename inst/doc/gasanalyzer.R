## ----setup, include = FALSE, eval = knitr::is_latex_output()------------------
#  knitr::opts_chunk$set(collapse = TRUE, comment = "#>", dev = "cairo_pdf")
#  # if you want to build vignettes, set NOT_CRAN to true in your .Renviron!
#  NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

## ----setup2, include = FALSE, eval = !knitr::is_latex_output()----------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
# if you want to build vignettes, set NOT_CRAN to true in your .Renviron!
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

## ----libs, message = FALSE----------------------------------------------------
library(gasanalyzer)
library(units)

## ----morelibs, message = FALSE, eval = NOT_CRAN, purl = NOT_CRAN--------------
library(ggplot2)
library(gridExtra)

## ----theming, include=FALSE, eval = NOT_CRAN, purl = NOT_CRAN-----------------
old_options <- options(digits = 2)

theme_pub <- theme_bw() + theme(
  plot.title = element_text(size = 14),
  plot.subtitle = element_text(size = 12),
  axis.text = element_text(size = 10),
  axis.ticks.x.top = element_blank(),
  axis.title = element_text(size = 12),
  legend.text = element_text(size = 11),
  legend.title = element_text(size = 12),
  panel.grid.minor = element_blank())
theme_set(theme_pub)

## ----o2data, message=FALSE----------------------------------------------------
exampleFiles <- system.file("extdata", package = "gasanalyzer")
# import calibration factors for the 6800 used here:
import_factory_cals(exampleFiles)
# load data
xlsxfile <- read_6800_xlsx(paste(exampleFiles, "lowo2.xlsx", sep = "/"))
txtfile <- read_6800_txt(paste(exampleFiles, "lowo2", sep = "/"))

## -----------------------------------------------------------------------------
# exclude list columns with equations, and divisions by zero
excludeCols <- vapply(xlsxfile, 
                      FUN=\(x) is.list(x) || 
                        all(is.infinite(x)) ||
                        all(is.nan(x)), TRUE)

# Interestingly, Leak.Fan reported in the txt file differs slightly
# from that in the xlsx. It is likely caused by rounding/averaging in the
# instrument firmware.
all.equal(txtfile[!excludeCols], xlsxfile[!excludeCols], tolerance = 1e-6)

## ----plot1, fig.height = 4, fig.width = 5.5, eval = NOT_CRAN, purl = NOT_CRAN----
# default points:
update_geom_defaults("point", list(shape = 1, size = 3))

# convenience function to change plot labels automatically:
nice_labels <- function(plt) {
  plt$labels <- var2label(plt$labels)
  plt
}

#generate some descriptive stats  and plotmath it with a call to var2label:
descr <- colMeans(xlsxfile[c("LeafQ.Qin", "Meas.Tleaf", "Const.Oxygen")])
dlist <- var2label(rep(names(descr), 2),
                   use_units = rep(c(FALSE, TRUE), each = 3),
                   val = c(rep(NA, 3), as.character(descr)),
                   group = c("", "")) |>
    setNames(letters[1:6]) |> unlist() |>
    substitute(a==d*","~b==e*","~c==f, env = _ ) |> as.expression()

AvsCi1 <- xlsxfile |> ggplot(aes(GasEx.Ci, GasEx.A)) + geom_point() +
  labs(title = dlist)

nice_labels(AvsCi1)

## -----------------------------------------------------------------------------
gswOld <- xlsxfile$GasEx.gsw
xlsxfile <- xlsxfile |> 
  transform(Const.K = 0.22) |>
  recalculate()
#max effect:
set_units(max(1 - as.numeric(xlsxfile$GasEx.gsw / gswOld)) * 100, "%")

## -----------------------------------------------------------------------------
#confirm that the slope of A vs ETR is not 1:
xlsxfile |> 
    lm(GasEx.A ~ I(FLR.ETR/4), data = _) |> coef() |> _[[2]]

## -----------------------------------------------------------------------------
# a linear regression:
lm1 <-
  xlsxfile |> transform(FLR.PhiQin_4 = FLR.phiPS2 * LeafQ.Qin / 4) |>
  lm(GasEx.A ~ FLR.PhiQin_4, data = _)

alphabeta <- lm1$coef[[2]]

## ----plot2, fig.height = 4, fig.width = 5.5, eval = NOT_CRAN, purl = NOT_CRAN----

AvsPhiQin1 <- ggplot(lm1$model, aes(x = .data[[names(lm1$model)[2]]],
                    y = .data[[names(lm1$model)[1]]])) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    labs(subtitle = bquote(alpha ~ beta == .(lm1$coef[[2]]) ~~~
                               R[L] == .(-lm1$coef[[1]])))

nice_labels(AvsPhiQin1)


## -----------------------------------------------------------------------------
# a long equation estimating the light absorption:
xlsxEqs <- xlsxfile$gasanalyzer.Equations[[1]]
print(xlsxEqs$LeafQ.alpha)

# add our custom alpha by dividing it through the used beta:
xlsxfile$Custom.alpha <- alphabeta / xlsxfile$Const.fPS2
# and change the equation to use the new value:
xlsxEqs <- modify_equations(xlsxEqs, LeafQ.alpha = \() {Custom.alpha})

# make a copy of the data and apply the new equations:
xlsxfileCopy <- xlsxfile |> 
  transform(gasanalyzer.Equations = list(xlsxEqs)) |>
  recalculate()

# confirm that the slope of A vs ETR is now 1:
xlsxfileCopy |> 
  lm(GasEx.A ~ I(FLR.ETR / 4), data = _) |> coef() |> _[[2]]


## ----tryeqs-------------------------------------------------------------------
RecalEqs <- create_equations(c("li6800", "default", "O2_correction"), 
                             LeafQ.alpha = \() {Custom.alpha})
xlsxfileRecal <- recalculate(xlsxfile, RecalEqs)

all.equal(xlsxfileRecal[names(xlsxfile)], xlsxfile, tol = 0.001)

## ----eqs----------------------------------------------------------------------
# providing a bad name will show a warning with valid names:
try <- create_equations("help")

# an equation set for gm, with a custom equation for resistance:
create_equations("gm_fluorescence", FLR.rm = \(x) {1 / FLR.gm})

## -----------------------------------------------------------------------------
xlsxfile21 <- xlsxfileRecal |> 
  transform(Const.Oxygen = set_units(21, "%")) |>
  recalculate(RecalEqs)

## ----plot3, fig.height = 3, fig.width = 6, eval = NOT_CRAN, purl = NOT_CRAN----
AvsCi21 <- xlsxfile21 |> ggplot() + 
  geom_point(aes(GasEx.Ci, GasEx.A, shape = "21% O2")) +
  geom_point(aes(xlsxfileRecal$GasEx.Ci,
                 xlsxfileRecal$GasEx.A, shape = "1% O2")) +
  scale_shape_manual(name = "Calculated at:",
                     values = c("1% O2" = 1, "21% O2" = 2)) +
  theme(legend.position = "inside", legend.position.inside = c(0.7, 0.3))

AvsPhiQin21 <-
  xlsxfile21 |>
  lm(GasEx.A ~ FLR.PhiQin_4, data = _) |> (\(x) {
    ggplot(x$model, aes(x = .data[[names(x$model)[2]]],
                        y = .data[[names(x$model)[1]]])) +
      geom_point() +
      stat_smooth(method = "lm", formula = y ~ x) +
      labs(subtitle = bquote(alpha ~ beta == .(x$coef[[2]]) ~~~
                            R[L] == .(-x$coef[[1]]))) +
      geom_point(data = xlsxfileRecal, shape = 2) 
    })()

grid.arrange(nice_labels(AvsCi21), nice_labels(AvsPhiQin21), ncol = 2)

## -----------------------------------------------------------------------------
xlsxfileRecal |> 
    lm(GasEx.A ~ I(FLR.ETR / 4), data = _) |> coef()

xlsxfile21 |> 
    lm(GasEx.A ~ I(FLR.ETR / 4), data = _) |> coef() 

## ----loadACi------------------------------------------------------------------
aciFiles <- list.files(exampleFiles, "aci\\d\\.csv", full.names = TRUE)

ACi <- lapply(aciFiles, \(x) {print(basename(x)); read_gfs(x, tz = "CEST")}) |>
   do.call("rbind", args = _)

data.frame(aggregate(list(RH = ACi$GasEx.RHcham), 
                     list(file = ACi$SysObs.Filename), mean),
           RH.SD = aggregate(list(ACi$GasEx.RHcham), 
                     list(ACi$SysObs.Filename), sd)[, 2]) |>
  knitr::kable()


## -----------------------------------------------------------------------------
gfsEqs <- create_equations(c("default", "gfs3000"))
ACi2 <- recalculate(ACi, gfsEqs)
all.equal(ACi$GasEx.Ci, ACi2$GasEx.Ci)

## -----------------------------------------------------------------------------
# we define some constants that specify the relative value of the variables
# on which we will do a sensitivity analysis:
SAset <- c("Const.calk", "Const.calg", "Const.calrhi", "Const.calt",
           "Const.calh")
# symbols for the plot header:
SAlabs <- c("italic(K)", "italic(g)[cw]", "RH[i]", 
            "italic(T)[leaf]", "'['*H[2]*O*']'[r]")
# set defaults to 1 (100%) and define sequences over which to permutate
ACi[SAset] <- 1
# a 5% range:
seq5 <- seq(-0.05, 0.05, length.out = 11)
seq5 <- seq5[seq5 != 0]

# make equations allowing to vary the variables by a relative value:
MSFeqs <- create_equations(c("default", "gfs3000", "cuticular_conductance"), 
                           Const.K = function() { Const.calk * Const.K },
                           Meas.Tleaf = function() { Const.calt * Meas.Tleaf },
                           Const.gcw = function() { Const.calg * Const.gcw },
                           Const.gcc = function() { Const.gcw / 20 },
                           Const.RHi = function() { Const.calrhi * Const.RHi },
                           Meas.H2Or = function() { Const.calh * Meas.H2Or },
                           Meas.H2Os = function() { Const.calh * Meas.H2Os })

# set gcw and use steady-state equations, then use the sequences
ACi[["Const.gcw"]] <- set_units(25, "mmol*m^-2*s^-1")
ACi[["Const.RHi"]] <- set_units(95, "%")
ACi[["SysConst.UseDynamic"]] <- FALSE

ACi.SA <- rbind(permutate(ACi, Const.calk = c(1, 1 + 20 * seq5)),
              permutate(ACi, Const.calg = 1 + 20 * seq5),
              permutate(ACi, Const.calrhi = 1 + seq5),
              permutate(ACi, Const.calt = 1 + seq5),
              permutate(ACi, Const.calh = 1 + seq5)) |>
  recalculate(MSFeqs)

## ----fitaci, eval = NOT_CRAN, purl = NOT_CRAN---------------------------------
library(photosynthesis)

aPars <- c("V_cmax", "J_max", "V_TPU", "R_d")

fitaci_per_group <- function(df, aPars, g1, ...) {
  grp <- c(g1, ...)
  outCols <- c(grp, aPars, "group", "Plots")
  splitOn <- paste("~", paste(grp, collapse = " + "))
  # drop category name and Filename for shorter labels:
  lbls <- vapply(strsplit(grp, ".", fixed = TRUE), \(x) rev(x), c("nm","gp")) |>
    paste0(": ") |>
    gsub("Filename: ", "", x = _, fixed = TRUE) 
  # layout:
  plotAdd <- list(do.call(labs, var2label(list(y = "GasEx.A", x = "GasEx.Ci"),
                                         TRUE)),
                  theme(legend.position = "none", plot.title.position = "plot",
                        plot.title = element_text(size = 8)))
  # translate names to what the photosynthesis package uses (sadly not essdive):
  out <- df[c("GasEx.A", "GasEx.Ci", "LeafQ.Qin", "Meas.Tleaf", 
              "Meas.Pa", "gasanalyzer.Permutate", grp)] |>
    setNames(c("A_net", "C_i", "PPFD", "T_leaf", "Pa", "group", grp)) |>
    drop_units() |> transform(T_leaf = T_leaf + 273.15) |>
    split(as.formula(splitOn), drop = TRUE) |>
    lapply(\(x) {
      lbl <- paste0(lbls, x[1, grp], collapse=", ")
      fit <- try(fit_aci_response(x, Pa = mean(x$Patm)))
      if (!inherits(fit, "try-error"))
        data.frame(x[1, c("group", grp)], 
                   fit[["Fitted Parameters"]],
                   Plots = I(list(fit$Plot + ggtitle(lbl) + plotAdd))
                   )[outCols]
      else 
        data.frame()
    }) 
  do.call(rbind, c(out, make.row.names = FALSE)) 
}

## ----eval = NOT_CRAN, purl = NOT_CRAN-----------------------------------------
# call the function to fit the curves:
CO2curves.SA <- ACi.SA |>
  fitaci_per_group(aPars, "SysObs.Filename", SAset)

# reference values are where the SAset variables equal 1:
refids <- rowSums(CO2curves.SA[SAset] == 1) == length(SAset)
# combine all relative values in one column based on group:
CO2curves.SA$gval <- unlist(CO2curves.SA[cbind(seq_len(nrow(CO2curves.SA)),
                                        match(CO2curves.SA$group,
                                              names(CO2curves.SA)))])
# value as relative change:
CO2curves.SA$gval <- 100 * (as.numeric(CO2curves.SA$gval) - 1)

# calc mean over the 3 reps:
CO2curves.SAm <- 
  lapply(names(CO2curves.SA[aPars]),
         \(x) {
           # relative change:
           rval <- 100 * as.numeric((CO2curves.SA[[x]] - 
                                       CO2curves.SA[[x]][refids]) /
                                      abs(CO2curves.SA[[x]][refids]))
           # grouping:
           byl <- as.list(CO2curves.SA[SAset])
           with(CO2curves.SA, 
                data.frame(nm = x, 
                           aggregate(cbind(gval), byl, head, n = 1),
                           group = aggregate(group, byl, head, n = 1)[["x"]],
                           val = aggregate(rval, byl, mean)[["x"]],
                           sd = aggregate(rval, byl, sd)[["x"]]))}) |>
  do.call(rbind, args = _)

## ----fig.height = 2.5, fig.width = 7, warning = FALSE, eval = NOT_CRAN, purl = NOT_CRAN----
grid.arrange(grobs = CO2curves.SA$Plots[refids], ncol = 3)

## ----eval=NOT_CRAN, fig.height = 4.5, fig.width = 7, purl=NOT_CRAN------------
# make lbls factors for good order and plotmath
CO2curves.SAm$group <- factor(CO2curves.SAm$group,
                              levels = SAset,
                              labels = SAlabs)
CO2curves.SAm$nm <- factor(CO2curves.SAm$nm,
                           levels = c("V_cmax", "J_max", "V_TPU", "R_d"),
                           labels = c("italic(V)[cmax]", "italic(J)[max]",
                                      "italic(V)[TPU]", "italic(R)[L]"))
# and plot:
CO2curves.SAm |>
  ggplot(aes(x = gval, y = val, group = nm)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotdash", color = "darkgrey") +
  geom_point() +
  geom_ribbon(aes(ymin = val - sd, ymax = val + sd), alpha = 0.1) +
  facet_grid(cols = vars(group), rows = vars(nm),
             scales="free", labeller = label_parsed) +
  xlab("Relative change in variable (%)") + 
  ylab("Relative change in parameter (%)") +
  theme(plot.title = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(14, "pt"))



## ----d13Crecalc---------------------------------------------------------------
# load data
isoFile <- file.path(exampleFiles, "d13C.tsv")
# read and recalculate
isotopes <- read_gasexchange(isoFile) |>
  recalculate(create_equations(c("default", "li6400")))

# vary f between 0 and 20, and use 2 sets of equations:
isotopes$d13CConst.e <- set_units(0, "permille")
isotopes$gasanalyzer.Equations <- list(NA)
isotopes <- isotopes |>
  permutate(d13CConst.f = seq(0, 20, 0.5)) |> 
  permutate(gasanalyzer.Equations = 
              c(connected = list(create_equations("d13C")),
                disconnected = list(create_equations(c("d13C", 
                                                       "d13C_dis"))))) |>
  recalculate() |>
  transform(factor.f = factor(d13CConst.f))

# point method gm is already calculated, average the results:
gmEst <- aggregate(list(gm = as.numeric(isotopes$d13C.gm)),
                   list(d13CConst.f = isotopes$d13CConst.f,
                        Model = isotopes$gasanalyzer.PermutateLabel),
                   FUN = function(x) 
                     c(d13C.gm = mean(x), d13C.gm.sem = sd(x)
                       / sqrt(length(x))))
gmEst <- do.call(data.frame, c(gmEst[1:2], unname(gmEst[3])))
gmEst$method <- "point"


## ----d13Cslopes, fig.height = 4, fig.width = 5.5, eval = NOT_CRAN, purl = NOT_CRAN----

d13Cslopes <- isotopes |>
  subset(gasanalyzer.PermutateLabel == "connected" &
           as.numeric(d13CConst.f) %in% c(7, 11, 16)) |>
  ggplot(aes(y = d13C.DeltaiDeltao, x = d13C.A_pCa, color = factor.f)) + 
  geom_point() +  
  geom_smooth(method = "lm", formula = y ~ x, aes(fill = factor.f)) +
    guides(fill = "none", 
         colour = guide_legend(override.aes = list(fill = NA))) +
  labs(color = var2label("d13CConst.f", TRUE)[[1]]) +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.10),
        legend.direction = "horizontal")

nice_labels(d13Cslopes)


## ----eval=NOT_CRAN, purl=NOT_CRAN---------------------------------------------
# The equation for the slope method is already in the data, so this function
# just evaluates that. Its not the fastest way, because recalculate adds
# overhead and unneccessary calcs, but simple:
fitfunc <- function(gm, ep) {
    # Di in the equation doesn't take gm into account, and assumes Ci=Cc
    # to get the true Delta we so we modify Ci to be Cc again:
    df$GasEx.Ci <-  df$GasEx.Ci - df$GasEx.A / 
      (set_units(gm, "mol*m^-2*s^-1*bar^-1") * df$Meas.Pa)
    # substitute estimated ep in data and recalc:
    df$d13C.ep <- set_units(ep, "permille")
    as.numeric(recalculate(df)$d13C.Deltai)
}

#this will take a while:
gmEstSlope <- isotopes |> 
  #we want to fit ep, so modify the equation list to not recalculate it:
  transform(gasanalyzer.Equations = 
              lapply(gasanalyzer.Equations,
                     \(x) modifyList(x, list(d13C.ep = NULL)))) |>
  #split by f and equation version, and call nls: 
  split(~ isotopes$d13CConst.f + isotopes$gasanalyzer.PermutateLabel) |>
  lapply(function(dat) {
           # data is moved to environment of fitfunc
           environment(fitfunc) <- list2env(list(df = dat))
           nlsfit <- nls(drop_units(d13C.Deltao) ~ fitfunc(gm, ep),
                         start = c(gm = 0.3, ep = -10), dat)
           cfs <- summary(nlsfit)$coefficients
           data.frame(d13CConst.f = dat$d13CConst.f[1],
                      Model = dat$gasanalyzer.PermutateLabel[1],
                      d13C.gm = cfs[1],
                      d13C.gm.sem = cfs[3], d13C.ep = cfs[2], 
                      d13C.ep.sem = cfs[4])}) |>
  do.call(rbind, args = _) |>
  transform(method = "slope")

gmEst <- rbind(gmEst, gmEstSlope[names(gmEstSlope) %in% names(gmEst)],
               make.row.names = FALSE)

#re-adds units:
units(gmEst$d13C.gm) <- deparse_unit(isotopes$d13C.gm)
units(gmEst$d13C.gm.sem) <- deparse_unit(isotopes$d13C.gm)


## ----fsens, fig.height = 4, fig.width = 7, eval = NOT_CRAN, purl = NOT_CRAN----

fsens <- ggplot(gmEst, aes(d13CConst.f, d13C.gm, 
                           ymin = d13C.gm - d13C.gm.sem,
                           ymax = d13C.gm + d13C.gm.sem, 
                           shape = method, color = Model)) + 
  geom_point(color = "black") +
  geom_ribbon(alpha = 0.1, color = NA) + guides(color = "none") +
  scale_shape_manual("Method:", values = c(1, 2)) + facet_wrap(vars(Model)) +
  theme(legend.position = "inside", legend.position.inside = c(0.12, 0.85),
        panel.grid.minor = element_blank()) 
nice_labels(fsens)

## ----include=FALSE, eval = NOT_CRAN, purl = NOT_CRAN--------------------------
options(old_options)

