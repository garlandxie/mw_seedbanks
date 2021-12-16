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
library(here)
library(googlesheets4)
library(visdat)
library(ggplot2)
library(dplyr)
library(forcats)
library(janitor)

# import ----

link <- "https://docs.google.com/spreadsheets/d/1O1Ll_PsW3qKwdZ_xnTrDKT_kGnpLGvtUKBQ73zvvBBM/edit?usp=sharing"
df_spr <- googlesheets4::read_sheet(link, sheet = "raw_data")

# check packaging ----

str(df_spr)
head(df_spr, n = 10)
tail(df_spr, n = 10)

# check missing data ----

visdat::vis_miss(df_spr)
visdat::vis_dat(df_spr)

# data clean -----

# re-assign species codes 
df_spr_tidy <- df_spr %>%
  janitor::clean_names() %>%
  mutate(spp_code = case_when(
    spp_code == "FRVE"  ~ "PONO", 
    spp_code == "ERCA"  ~ "COCA", 
    spp_code == "PACA"  ~ "PAMI", 
    spp_code == "AMSP"  ~ "VETH", 
    spp_code == "CAREX" ~ "ANGE", 
    spp_code == "SEVI"  ~ "ANGE",
    TRUE ~ spp_code)
  )
  
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
