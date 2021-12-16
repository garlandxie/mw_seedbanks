################################################################################
# Accompanying code for the followng research project: 
#  Seedbanks in the Meadoway (Scarborough, Canada)
#
#
#  Corresponding authors for this script:  
#   Garland Xie      (1)
#
# Affiliations: 
#   (1) Department of Biological Sciences, 
#       University of Toronto Scarborough,
#       1265 Military Trail, Toronto, ON, M1C 1A4, Canada
#       email: garland.xie@mail.utoronto.ca, 
#              nicholas.sookhan@mail.utoronto.ca
#              scott.macivor@mail.utoronto.ca
#
# Purpose of this R script: to process taxonomic data for the seedling emergents
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# libraries ----
library(here)          # for creating relative file-paths
library(googlesheets4) # for parsing google spreadsheet files
library(visdat)        # for visualizing missing data 
library(dplyr)         # for manipulating data 
library(janitor)       # for cleaning colum names in a machine-readable format 

# import ----

link <- "https://docs.google.com/spreadsheets/d/1r-T9lY1Osjez8SKKoCUHRi5ftbHBvgFn2Xk-IPqc0z4/edit?usp=sharing"
taxon_df <- googlesheets4::read_sheet(link, sheet = "raw_data")

# check packaging ---

glimpse(taxon_df)
head(taxon_df, n = 5)
tail(taxon_df, n = 5)

# check for taxonomic synyonyms ----

# use taxize R package

# save to disk -----

write.csv(
  x = taxon_df,
  here("data", "analysis_data", "taxonomic_info.csv")
)

