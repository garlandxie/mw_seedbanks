################################################################################
# Accompanying code for the followng research project: 
#   Seedbanks in the Meadoway (Scarborough, Canada)
#
#
# Corresponding authors for this script:  
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
# Purpose of this R script: to process data for the spring seedbank season
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# libraries ----
library(here)          # for creating relative file-paths
library(visdat)        # for visualizing missing data 
library(ggplot2)       # for visualizing data 
library(dplyr)         # for manipulating data 
library(janitor)       # for cleaning colum names in a machine-readable format 

# import ----

df_spr <- read.csv(here("data", "input_data", "seed_bank_spring.csv"))
  
# check packaging ----

str(df_spr)
head(df_spr, n = 10)
tail(df_spr, n = 10)

# check missing data ----

visdat::vis_miss(df_spr)
visdat::vis_dat(df_spr)

# data clean -----

# summarize
df_spr_summ <- df_spr_tidy %>%
  group_by(season, section, site, treatment, plot, spp_code) %>%
  summarize(total_abund = sum(abund, na.rm = TRUE)) %>%
  ungroup()

# sanity checks
glimpse(df_spr_summ)
vis_miss(df_spr_summ)

# write to disk ----

write.csv(
  x = df_spr_summ, 
  file = here("data", "analysis_data", "spring_seedbank.csv")
)
