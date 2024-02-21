# Seedbanks in the Meadoway (Scarborough, Canada)

The following documentation of this README file follows the TIER 4.0 Protocol (https://www.projecttier.org/tier-protocol/protocol-4-0/root/)

## Software and version

Code to reproducible the analysis and figures is written in the R programming language (version 4.3.2). 
The developer (Garland Xie) typically writes and runs the code using macOS Big Sur 11.6

## Abstract 

Designed urban meadows can rehabilitate wildlife habitat and increase floral biodiversity in formerly degraded spaces such as beneath electrical power lines where trees are not permitted, and mowing is frequent. Monitoring seed banks can evaluate the efficacy of restoration efforts, including proliferation of native seed mixes or the reduction in non-curated seeds from nearby plant populations. Using soil seed banks, we investigate how above-ground restoration management throughout a 16-km infrastructure corridor in Toronto, Canada impacts seedling recruitment of native seed mix, spontaneous. and non-native plant species. Soil cores were taken in spring and fall from 80 plots across nine meadow locations at two different restoration stages: 1) newly-established sites that received a seed bank preparation of rototilling, cover crops, and application of native seeds through seed drilling, and 2) restored (5-8 years) sites that alternates between two maintenance treatments (unmown and mown). Soil samples were transplanted into greenhouse trays and surveyed for 100 days. A total of 11,373 seedlings representing 93 species of forbs, grasses, and woody plants were surveyed. Restored meadow habitats had at least 30% more species, 3 times more seeds, and the proportion of native seed mix species was 9 times higher compared to newly-established habitats. Spring sampling had 0.8 times more seedlings and 1.5 times higher proportion of invasive species in restored plots compared to fall. Our findings show the practical implications of seed bank dynamics in designed urban meadows, where the optimal timing to evaluate potential recruitment of native and invasive seedlings is in spring.

## Folder structure 

- Data
  - [InputData](data/input_data): all raw data 
    - [Metadata](data/input_data/metadata): containing information about the sources and contents of InputData files
  - [AnalysisData](data/analysis_data): data that is in a format suitable for the analysis   
  - [IntermediateData](data/intermediate_data): partially processed data to use in subsequent data analyses
    
- Scripts
  - [ProcessingScripts](scripts/processing_scripts): The commands in these scripts transform Input Data Files into Analysis Data Files
  - [AnalysisScripts](scripts/analysis_scripts): The commands in these scripts generate the results in the manuscript
- Output
  -  [DataAppendixOutput](output/data_appendix_output): contains files with figures, tables and other output presented in the Data Appendix
  -  [Results](output/results): contains files with figures, tables, and other output presented in the manuscript

## Instructions for reproducing the results

To use the code in this repository to reproduce the manuscript's results,
please follow the following steps:
1. `git clone` this repository or download it as a zip folder
2. Open `Rstudio`, go to `file > Open project` and open the `Msc.Rproj`
Rproject associated with this repository
3. Run `renv::restore()` in your R console. Requires `renv` package (see [THIS](https://rstudio.github.io/renv/articles/renv.html) vignette)

