# gasanalyzer 0.4.3
* Minor update
  - Work around undefined behaviors when using empty vectors with units
  - Add methods for boundary layer conductance calculations
  - Add missing detailed help for create_equations
  - Fix another incorrect timing related potential test failure
  
# gasanalyzer 0.4.2
* Minor update
  - Make xlsx import handle spreadsheet errors better
  - Fix an bug with equation list when UseFlags were not the last list entry
  - Fix CRAN check error with a test
  
# gasanalyzer 0.4.1
* Minor update
  - setup test infrastructure (not many useful unit tests added yet)
  - bug fixes to gas-exchange equations (especially related to cuticular conductance)
  - Add option to modulate the relative humidity in the leaf
  - reorder some variables in more logical categories 
  - Fix loading of some LI6800 files
  - Make file loading in examples more platform independent
  - Updated description and refer to the publication about the package
  
# gasanalyzer 0.4.0

* Minor update
  - Add support for CIRAS-4 (PP Systems) instruments
  - Many small code style fixes
  - Minor corrections to how stomatal ratios are handled in MSF equations
  - Corrected GFS3000 boundary conductance equation
  - Additional checks are done now when correcting for O2 

# gasanalyzer 0.3.4

* CRAN resubmission:
  - Add return value even if no return value
  - Avoid making code comments, dontrun and interactive in examples.
  - Take care to reset user options
  - Added some references to the description

# gasanalyzer 0.3.3

* Initial CRAN submission.
