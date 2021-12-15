# mw_seedbanks
Perennial seedbank surveys along a 16km stretch of urban hydrocorridors 

The following documentation of this README file follows the TIER 4.0 Protocol (https://www.projecttier.org/tier-protocol/protocol-4-0/root/)

### Using the code in this repository:
Code to reproducible the analysis and figures is written in the R programming language. 

To use the code in this repository to reproduce the manuscript's results,
please follow the following steps:
1. `git clone` this repository or download it as a zip folder
2. Open `Rstudio`, go to `file > Open project` and open the `Msc.Rproj`
Rproject associated with this repository
3. Run `renv::restore()` in your R console. Requires `renv` package (see [THIS](https://rstudio.github.io/renv/articles/renv.html) vignette)

## Folder structure 

- Report
- Data
  - InputData: all raw data 
    - Metadata: containing information about the sources and contents of InputData files
  - AnalysisData: data that is in a format suitable for the analysis   
  - IntermediateData: partially processed data to use in subsequent data analyses
- Scripts
  - ProcessingScripts: The commands in these scripts transform Input Data Files into Analysis Data Files
  - DataAppendixScripts: The commands in theses scripts produce the figures, tables, and descriptive statistics in the Data Appendix
  - AnalysisScripts: The commands in these scripts generate the results in the Rmarkdown report
- Document: contains the master script that reproduces the Results of your project by executing all the other scripts, in the correct order
- Output
  -  DataAppendixOutput: contains files with figures, tables and other output presented in the Data Appendix
  -  Contains files with figures, tables, and other output presented in your report

