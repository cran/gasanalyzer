---
title: "gasanalyzer"
output: html_document
---
<!-- badges: start -->

<!-- badges: end -->

Provides support for reading and preprocessing data from various portable photosynthesis systems.
Although there are other tools for importing data from LI-COR instruments[^1]<sup>,</sup>[^2]<sup>,</sup>[^3]<sup>,</sup>[^4],
they are more limited in scope. Gasanalyzer aims to import data from a wider range of instruments and formats, convert the data
to a tidy format, and also support the calculation of derived quantities from such data.

To achieve the latter goal, the package allows importing data and equations stored in xlsx files, provides sets of pre-defined equations to
calculate physiologically relevant derived variables, and allows users to add additional equations for preprocessing and analysis

### Installation
Gasanalyzer can be found on [CRAN](https://cran.r-project.org/package=gasanalyzer), which greatly simplifies installation
(if you are using an up-to-date R install, it can for example be found in RStudio under Tools, Install Packages).

It is also possible to install the development version of the package from gitlab. Please make sure to run the following command
to install all requirements first:

```r
install.packages(c("units", "stringi", "jsonify", "xml2", "tidyxl", "tibble", "vctrs", "devtools"))
```

Then, download the zip (or tar.gz, etc) from gitlab and run:

```r
devtools::install_local("gasanalyzer-master.zip")
```
You may need to adopt the command above to point to the correct folder and name of the file you just downloaded.
R may give a warning about missing Rtools, however, Rtools is currently not needed to install this package.

### Documentation

Portable photosynthesis systems are infra-red gas-analyzers for
measuring gas-exchange characteristics on plant leaves. They typically
measure $CO_{2}$ and $H_{2}O$ mol fractions, gas flow and various other
relevant parameters (temperature, light intensity, fan speed, pressure).
These measurements are combined with user-defined parameters (leaf area
or weight, stomatal ratios, oxygen concentration) and used by the
instrument firmware or external software to calculate physiological
relevant traits such as the rate of photosynthesis, evapotranspiration,
or intercellular $CO_{2}$ concentrations. These calculations have been described
in scientific publications[^5]<sup>,</sup>[^6] and user manuals for the
instruments[^7]<sup>,</sup>[^8]<sup>,</sup>[^9]. Unfortunately, there are some
differences in the approaches, assumptions and terminology used by the different
instruments. This makes it difficult to compare the results.

*A unified way of dealing with gas-exchange data would benefit the research in this field.*
Moreover, a change in assumptions, configuration or externally measured data (leaf area,
stomatal ratio, $O_2$) makes it often necessary to recalculate gas-exchange data. Although some
vendor-provided spreadsheets provide some options to recalculate data, they are
not available for all instruments, and limited in scope and usability. Thus, *A tool to reliably
recalculate gas-exchange data is currently not available*. To this end, this *R*
package uses a unified set of symbols and equations for recalculating gas-exchange data.

An advantage of using *R* is that it allows us to read and modify many datafiles
as part of a single scripted procedure. Afterwards the data can be quickly
summarized, analyzed and plotted. It ensures a repeatable and traceable pipeline
to turn raw measurements into analyzed results.

Please see the package vignette for more detailed information and examples on how to use this package.

[^1]: https://github.com/erikerhardt/RLicor
[^2]: https://github.com/poales/readLicorData
[^3]: https://github.com/muir-lab/licorer
[^4]: https://github.com/PaulESantos/licor6400r
[^5]: von Caemmerer and Farquhar (1981). Some relationships between the
    biochemistry of photosynthesis and the gas exchange of leaves.
    Planta 153, 376--387. doi: 10.1007/bf00384257
[^6]: Márquez, Stuart-Williams and Farquhar. An improved theory for
    calculating leaf gas exchange more precisely accounting for small
    fluxes. Nat. Plants 7, 317--326 (2021). doi:
    10.1038/s41477-021-00861-w
[^7]: <https://www.licor.com/env/support/LI-6400/topics/system-description.html#Equation>
[^8]: <https://www.licor.com/env/support/LI-6800/topics/equation-summary.html>
[^9]: <https://www.walz.com/files/downloads/manuals/gfs-3000/GFS-3000_Manual_9.pdf>
[^10]: <https://ppsystems.com/download/technical_manuals/80150-1-CIRAS-4_Operation_V1.3.pdf>


```
