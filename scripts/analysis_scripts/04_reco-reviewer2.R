################################################################################
# Accompanying code for the following research project: 
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
#
#
# Purpose of this R script: to reanalyze the species composition data
# based on reviewer concerns from the Restoration Ecology submission

# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html


# re-run the original analysis for species composition, but split 
# it into two separate constrained correspondence analyses:

# (1) constrain restoration stage, but conditioning on sites + sampling season

# (2) remove newly-established sites (to acknowledge the artefact), 
# and constrain sampling season, while conditioning on site + restoration stage

# libraries --------------------------------------------------------------------
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data
library(ggplot2)    # for visualizing data
library(vegan)      # for doing community ecology analyses
library(ggrepel)    # for repelling text in ggplot2 figures
library(rdacca.hp)  # for calculating variance partitioning 

# import -----------------------------------------------------------------------

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

# clean data -------------------------------------------------------------------

sb <- rbind(fall_sb, spr_sb)

comm_matrix <- sb %>%
  group_by(season, treatment, site_code) %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = total_abund) %>%
  ungroup() %>%
  mutate(across(.cols = AMRE:LOPE, ~ tidyr::replace_na(.x, 0))) %>%
  mutate(
    season = factor(season, levels = c("Spring", "Fall")),
    treatment = factor(treatment, levels = c("TIL", "RES", "MOW"))
    )

# run RDA (with newly-established sites) ---------------------------------------

## prep ------------------------------------------------------------------------

spp_comp_w_til <- comm_matrix %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

env_w_til <- comm_matrix %>%
  dplyr::select(season, site_code, treatment)

## run the model ---------------------------------------------------------------

rda_w_til <- rda(
  spp_comp_w_til ~ treatment*season + site_code, 
  data = comm_matrix
  )

## tests of significance -------------------------------------------------------

anova.cca(rda_w_til, permutations = how(nperm = 999), by = "margin")
anova.cca(rda_w_til, permutations = how(nperm = 999), by = "axis")

## variance partitioning --------------------------------------------------------

vp_w_til <- rdacca.hp(
  dv = spp_comp, 
  iv = env_df, 
  method = "RDA",
  type = "adjR2"
)

# run RDA (without newly-established sites) ------------------------------------

## prep ------------------------------------------------------------------------

spp_comp_no_til <- comm_matrix %>%
  filter(!treatment %in% c("TIL")) %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

env_no_til <- comm_matrix %>%
  filter(!treatment %in% c("TIL")) %>%
  dplyr::select(season, site_code, treatment)

## run the model ---------------------------------------------------------------

rda_no_til <- rda(Y ~ treatment*season + site_code, data = comm_matrix)

## tests of significance -------------------------------------------------------

anova.cca(rda_no_til, permutations = how(nperm = 999), by = "margin")
anova.cca(rda_no_til, permutations = how(nperm = 999), by = "axis")

## variance partitioning --------------------------------------------------------

vp_no_til <- rdacca.hp(
  dv = spp_comp_no_til, 
  iv = env_df, 
  method = "RDA",
  type = "adjR2"
)

