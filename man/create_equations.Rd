% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equations.R
\name{create_equations}
\alias{create_equations}
\title{Create a list of equations for recalculating gasanalyzer data.}
\usage{
create_equations(useflags = "default", ...)
}
\arguments{
\item{useflags}{character vector with the type of equations to create
(see Details). Leave empty to obtain the default set. An unknown flag
returns an empty list, and a warning listing all valid flags.}

\item{...}{custom equations. the arguments must tagged function expressions.
Note that the function body must be wrapped in curly brackets. The tags
will be matched against the names of a data frame when applying the return
value with \code{\link[=recalculate]{recalculate()}}.}
}
\value{
A list of language objects with equations
}
\description{
This function creates a list of equations that can be used to recalculate
gas-exchange data by passing the resulting object to the \code{\link[=recalculate]{recalculate()}}
method. Various \code{useflags} can be defined to tune the equations.
In addition, custom equations can be defined as arguments. Note that
the calculations may fail if commons are missing in the gas-exchange data.
}
\details{
The \code{useflags} argument currently supports several pre-defined sets
of equations.
\subsection{Gas-Exchange and fluorescence models}{

\describe{
\item{\code{default}}{
A set with the most commonly-used derived quantities in the gas-exchange
(\code{GasEx}) and chlorophyll fluorescence (\code{FLR}) categories. The
equations are described in detail by von Caemmerer and Farquhar (1981)
and LI-COR Biosciences Inc (2022).
}
\item{\code{GoffGratch1946}}{
Replaces or adds equations related to the calculation of
the saturated water pressure of the leaf (\code{GasEx.SVPleaf}) and chamber
air (\code{GasEx.SVPcham}). The calculation is based on that described by
Goff and Gratch (1946) and only takes the temperature into account.
}
\item{\code{Buck1981}}{
Replaces or adds equations for (\code{GasEx.SVPleaf}) and (\code{GasEx.SVPcham}).
based on the description by Buck (1981). Takes temperature and
pressure into account.
}
\item{\code{Buck1996}}{
Replaces or adds equations for (\code{GasEx.SVPleaf}) and (\code{GasEx.SVPcham})
based on the description by Buck (1996). Takes temperature and
pressure into account.
}
\item{\code{cuticular_conductance}}{
Replaces the equations related to CO2 and H2O conductance
and substomatal CO2 (\code{GasEx.gtw}; \code{GasEx.gsw}; \code{GasEx.gtc}; \code{GasEx.Ci})
with versions that take into account cuticular conductance
(Márquez et al., 2021, 2023). Requires manually specifying the
cuticular conductance to water (\code{Const.gcw}) and CO2 (\code{Const.gcc}).
}
\item{\code{boundary_conducance}}{
Replaces or adds equations related to the boundary
layer conductance (\code{GasEx.gbw}) and leaf temperature derived from
the energy balance equations (\code{GasEx.TleafEB}) with a version that
is valid for a sample with no stomata. This is typically used when
estimating \code{GasEx.gbw} using filter paper (Parkinson, 1985). Note
that the leaf thermocouple (\code{Meas.Tleaf}) is used to estimate air
temperature (and should therefore not touch the sample). This can be
overridden by a custom equation. The sample emissivity (\code{Const.eps}),
the ratio between heat and water conductance for the chamber air
(\code{Const.ra_rv}) and the fraction of the area exchanging radiative and
sensible heat (\code{Const.asH}) can be adjusted. For solving the energy
balance equations, the method described by Bristow (1987) is used.
}
\item{\code{gm_fluorescence}}{
Adds derived variables for mesophyll conductance
(\code{FLR.gm}) and chloroplast CO2 mole fractions (\code{FLR.Cc}) based
on gas exchange and chlorophyll fluorescence (Harley et al., 1992). It is
strongly recommended to first calibrate the electron transport rates
reported in the \code{FLR.ETR} column. In addition, the equations require
adding columns for the respiration rate in the light (\code{Const.RL}),
and the CO2 photo-compensation point (\code{Const.GammaStar}).
}
}
}

\subsection{Isotope models}{

\describe{
\item{\code{d13C}}{
Adds derived variables related to the stable carbon isotope
discrimination model for C3 plants (Farquhar and Cernusak, 2012;
Evans and von Caemmerer, 2013) aimed at estimating the mesophyll
conductance \code{d13C.gm} and chloroplast CO2 mole fraction
(\code{d13C.Cc}). Requires additional columns with data on the carbon
isotope composition in sample and reference air
(\code{d13CMeas.delta13CO2s}; \code{d13CMeas.delta13CO2r}) and in the air
where the plants were grown (\code{d13CConst.delta13CO2g}), and values
for respiration in the light (\code{Const.RL}) and the CO2
photo-compensation point (\code{Const.GammaStar}).
}
\item{\code{d13C_dis}}{
Requires the \code{d13C} flag, but modifies the
modeled carbon isotope discrimination \code{d13C.Deltai} and
\code{d13C.gm} using the assumption that the carbon pools of respiration
and assimilation are disconnected (as described by Busch et al., 2020).
}
\item{\code{d13C_e_Busch2020}}{
Requires the \code{d13C} flag, but modifies the calculation
of the effective respiratory fractionation (\code{d13C.ep}) to
better take into account the effect of the growth conditions
(Busch et al. 2020). Additionally requires a value for the observed
discrimination against 13CO2 under growth conditions
(\code{d13CConst.Deltag}).
}
}
}

\subsection{Instrument related calculations}{

These implement calculations as close as possible to those used in the
instrument firmware and are typically related to the specific configurations
or instrument design. In addition, some methods are provided for
recalculating low-level instrument variables.
\describe{
\item{\code{li6400}}{
Definitions for boundary layer conductance, temperature, and light
related variables (\code{GasEx.gbw}; \code{GasEx.Rabs};
\code{GasEx.TairCnd}; \code{LeafQ.alpha}) specific to  LI-6400 / LI-6400XT
instruments (LI-COR Biosciences Inc, 2011).
}
\item{\code{li6800}}{
Definitions for boundary layer conductance, light, leakage,
temperature and light related variables (\code{GasEx.gbw};
\code{GasEx.Rabs}; \code{GasEx.TairCnd}; \code{LeafQ.alpha};
\code{LeafQ.Qin}; \code{LeafQ.Conv}; \code{Leak.Fan};
\code{Leak.CorrFact}; \code{GasEx.Ca}) specific to the LI-6800 instrument
(LI-COR Biosciences Inc, 2022).
}
\item{\code{ciras4}}{
Light absorptance (\code{LeafQ.alpha}) and the conversion between
photon flux density and irradiance (\code{LeafQ.Conv}) is calculated by
taking into account the effect of the different light sources that can be
used with the CIRAS-4 instrument (PP Systems, 2024).
}
\item{\code{gfs3000}}{
Boundary layer conductance (\code{GasEx.gbw}) and light sensor
(\code{LeafQ.Qin}) calculations adjusted for the default chamber of the
GFS-3000 instrument (Heinz Walz GmbH, 2019).
}
\item{\code{gfs3000_light_bot}}{
Requires gfs3000, but modifies \code{LeafQ.Qin} to indicate that the
bottom light sensor of the default chamber was used to quantify the
light intensity incident on the leaf (Heinz Walz GmbH, 2019).
}
\item{\code{match}}{
Takes a previously stored offset between sample and reference
analyzers into account when recalculating water and CO2 mole fractions
(\code{Meas.H2Os}; \code{Meas.CO2s}) as described by LI-COR Biosciences
Inc (2022), Heinz Walz GmbH (2019) and PP Systems (2024). The corrected
mole fractions are already stored in the data and therefore this
calculation is usually not required. However, these equations are needed
when recalculating mole fractions from lower level data (see the
\code{raw} and \code{O2_correction} flags).
}
\item{\code{raw}}{
Recalculates CO2 and H2O mole fractions from such low-level variables
(\code{Meas.CO2r}; \code{Meas.H2Or}; \code{Meas.CO2a}; \code{Meas.H2Oa})
as described by LI-COR Biosciences Inc (2022). Currently only
implemented for the LI-6800 because low-level instrument data are
required. Requires storing raw data and the availability of factory
calibration files. Requires and enables the \code{match} set.
}
\item{\code{O2_correction}}{
Recalculates CO2 and H2O mole fractions (\code{Meas.CO2r};
\code{Meas.H2Or};  \code{Meas.CO2a};  \code{Meas.H2Oa}) at a potentially
different oxygen concentration ( \code{Const.Oxygen}). Currently only
implemented for the LI-6800 (LI-COR Biosciences Inc, 2022) and GFS-3000
(K. Siebke, Heinz Walz GmbH, personal communication). For the LI-6800,
this requires loading of factory calibration files
(\code{\link[=import_factory_cals]{import_factory_cals()}}). Requires and automatically enables the
\code{match} flag.
}
}
}
}
\examples{
exampledir <- system.file("extdata", package = "gasanalyzer")

# import factory calibration for example data:
import_factory_cals(exampledir)

# read data from a txt file:
li6800 <- read_6800_txt(file.path(exampledir, "lowo2"))

# passing an invalid flags shows which flags are valid:
\donttest{create_equations("help")}

# create a default set of gas-exchange equations, for the 6800, but overwrite
# the default calculation of leaf light absorption with a custom value:
Eqs <- create_equations(c("default", "li6800"), LeafQ.alpha = \() {0.86})

#apply:
li6800_recalc <- recalculate(li6800, Eqs)

li6800$LeafQ.alpha
li6800_recalc$LeafQ.alpha
}
\references{
\itemize{
\item Bristow KL. 1987. On solving the surface energy balance equation for
surface temperature. Agricultural and Forest Meteorology \strong{39}:49–54.
\item Buck AL. 1981. New equations for computing vapor pressure and enhancement
factor. Journal of Applied Meteorology \strong{20}:1527–1532.
\item Buck AL. 1996. Buck research CR-1A user’s manual. Boulder, CO: Buck
Research Instruments.
\item Busch FA, Holloway-Phillips M, Stuart-Williams H, Farquhar GD. 2020.
Revisiting carbon isotope discrimination in C3 plants shows respiration
rules when photosynthesis is low. Nature Plants \strong{6}:245–258.
\url{https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf}
\item von Caemmerer S, Farquhar GD. 1981. Some relationships between the
biochemistry of photosynthesis and the gas exchange of leaves.
Planta \strong{153}:376–387.
\item Evans JR, von Caemmerer S. 2013. Temperature response of carbon isotope
discrimination and mesophyll conductance in tobacco. Plant, Cell &
Environment \strong{36}:745–756.
\item Farquhar GD, Cernusak LA. 2012. Ternary effects on the gas exchange of
isotopologues of carbon dioxide. Plant, Cell & Environment
\strong{35}:1221–1231.
\item Goff JA, Gratch S. 1946. Low-pressure properties of water from -160° F to
212° F. Transactions of the American Society of Heating and Ventilating
Engineers \strong{52}:95–121.
\item Harley PC, Loreto F, Di Marco G, Sharkey TD. 1992. Theoretical
considerations when estimating the mesophyll conductance to CO2 flux by
analysis of the response of photosynthesis to CO2. Plant Physiology
\strong{98}:1429–1436.
\item Heinz Walz GmbH. 2019. Portable gas exchange fluorescence system GFS-3000.
Handbook of operation. 9th edition.
\url{https://www.walz.com/files/downloads/manuals/gfs-3000/GFS-3000_Manual_9.pdf}
\item LI-COR Biosciences Inc. 2011. Using the LI-6400/LI-6400XT portable
photosynthesis system. Version 6.2.
\url{https://www.licor.com/env/support/LI-6400/manuals.html}
\item LI-COR Biosciences Inc. 2022. Using the LI-6800 portable photosynthesis
system. Version 2.1.
\url{https://www.licor.com/env/support/LI-6800/manuals.html}
\item Márquez DA, Stuart-Williams H, Farquhar GD. 2021. An improved theory for
calculating leaf gas exchange more precisely accounting for small fluxes.
Nature Plants \strong{7}:317–326.
\item Márquez DA, Stuart-Williams H, Cernusak LA, Farquhar GD. 2023. Assessing the
CO2 concentration at the surface of photosynthetic mesophyll cells.
New Phytologist \strong{238}:1446–1460.
\item Parkinson KJ. 1985. A simple method for determining the boundary layer
resistance in leaf cuvettes. Plant, Cell & Environment \strong{8}: 223–226.
\item PP Systems. 2024. CIRAS-4 portable photosynthesis system. Operation manual.
Amesbury, MA: PP Systems. Version 1.3.
}
}
\seealso{
\code{\link[=read_6800_equations]{read_6800_equations()}}
}
