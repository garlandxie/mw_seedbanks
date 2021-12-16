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
# Purpose of this R script: to process data for the fall seedbank season
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# libraries ----
library(here)
library(googlesheets4)
library(visdat)
library(ggplot2)
library(dplyr)
library(forcats)
library(janitor)

# import ----

link <- "https://docs.google.com/spreadsheets/d/1SWlk5eWdk3IOMFS9nv61p1HW-Tw4K3yR-tAIkAMlzs4/edit?usp=sharing"
df_fall <- googlesheets4::read_sheet(link, sheet = "raw_data")

# check packaging ----

str(df_fall)
head(df_fall, n = 10)
tail(df_fall, n = 10)

# check missing data ----

visdat::vis_miss(df_fall)
visdat::vis_dat(df_fall)

# data cleaning -----

# summarize
df_fall_summ <- df_sb %>%
  janitor::clean_names() %>%
  group_by(season, section, site, treatment, plot, spp_code) %>%
  summarize(total_abund = sum(abund, na.rm = TRUE)) %>%
  ungroup()

# sanity checks
glimpse(df_fall_summ)
vis_miss(df_fall_summ)

# write to disk ----

write.csv(
  x = df_fall_summ, 
  file = here("data", "analysis_data", "fall_seedbank.csv")
)