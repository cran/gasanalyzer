---
title: "gasanalyzer"
output: html_document
---
<!-- badges: start -->

<!-- badges: end -->

Provides support for reading and preprocessing data from various portable
photosynthesis systems. Includes options to parse xlsx equations from the data
files, for recalculating data based on sets of pre-defined equations, or based on
user-specified equations.

### Installation
GasAnalyzeR can be found on [CRAN](https://cran.r-project.org/package=gasanalyzer), which greatly simplifies installation (if you are using an up-to-date R install, it can for example be found in RStudio under Tools, Install Packages).

It is also possible to install the development version of the package from gitlab. Please make sure to run the following command to install all requirements first:

```r
install.packages(c("units", "stringi", "jsonify", "xml2", "tidyxl", "tibble", "vctrs", "devtools"))
```

Then, download the zip (or tar.gz, etc) from gitlab and run:

```r
devtools::install_local("gasanalyzer-master.zip")
```
You may need to adopt the command above to point to the correct folder and name of the file you just downloaded. R may give a warning about missing Rtools , however, Rtools is currently not needed to install this package.

### Documentation

Portable photosynthesis systems are infra-red gas-analyzers for
measuring gas-exchange characteristics on plant leaves. They typically
measure $CO_{2}$ and $H_{2}O$ mol fractions, gas flow and various other
relevant parameters (temperature, light intensity, fan speed, pressure).
These measurements are combined with user-defined parameters (leaf area
or weight, stomatal ratios, oxygen concentration) and used by the
instrument firmware or external software to calculate physiological
relevant traits such as the rate of photosynthesis, evapotranspiration,
or intercellular $CO_{2}$ concentrations.

These calculations are described in scientific publications[^1]<sup>,</sup>[^2]
and user manuals for the instruments[^3]<sup>,</sup>[^4]<sup>,</sup>[^5]. 
Some instruments optionally save data in the form of spreadsheet documents that 
contain the used equations. Unfortunately, these sources often made different 
assumptions for deriving gas-exchange equations. Moreover, they all use very 
different terminology and symbols. This makes it difficult to compare the 
different approaches. *A unified way of dealing with gas-exchange data would 
benefit the research in this field.* Moreover, a change in assumptions, 
configuration or externally measured data (leaf area, stomatal ratio, $O_2$)
makes it necessary to recalculate gas-exchange data. Although the 
vendor-provided spreadsheets provide some options to recalculate data, they are
not available and limited in scope and usability. *A tool to reliably 
recalculate gas-exchange data is currently not available*. To this end, this 
package uses a unified set of symbols and equations for gas-exchange data. 

An advantage of using *R* is that it allows us to read and modify many datafiles
as part of a single scripted procedure. Afterwards the data can be quickly
summarized, analyzed and plotted. It ensures a repeatable and traceable pipeline
to turn raw measurements into analyzed results.

Please see the package vignette for more detailed information.

[^1]: von Caemmerer and Farquhar (1981). Some relationships between the
    biochemistry of photosynthesis and the gas exchange of leaves.
    Planta 153, 376--387. doi: 10.1007/bf00384257
[^2]: Márquez, Stuart-Williams and Farquhar. An improved theory for
    calculating leaf gas exchange more precisely accounting for small
    fluxes. Nat. Plants 7, 317--326 (2021). doi:
    10.1038/s41477-021-00861-w
[^3]: <https://www.licor.com/env/support/LI-6400/topics/system-description.html#Equation>
[^4]: <https://www.licor.com/env/support/LI-6800/topics/equation-summary.html>
[^5]: <https://www.walz.com/files/downloads/manuals/gfs-3000/GFS-3000_Manual_9.pdf>
[^6]: <https://www.licor.com/env/support/LI-6800/topics/symbols.html>

```
