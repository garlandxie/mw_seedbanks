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

# run partial RDA to constrain restoration stage -------------------------------

## prep ------------------------------------------------------------------------

Y <- comm_matrix %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

## run the model ---------------------------------------------------------------

rda_trt_sn <- rda(Y ~ treatment*season + site_code, data = comm_matrix)
rda_sn <- rda(Y ~ season + site_code, data = comm_matrix)
rda_trt <- rda(Y ~ treatment + site_code, data = comm_matrix)
rda_int <- rda(Y ~ treatment*season, data = comm_matrix)

# variance partitioning --------------------------------------------------------

# get adjusted R-squared values for all RDA models
fra_d_R2 <- RsquareAdj(rda_trt_sn)$adj.r.squared
fra_b_R2 <- RsquareAdj(rda_sn)$adj.r.squared
fra_a_R2 <- RsquareAdj(rda_trt)$adj.r.squared
fra_c_R2 <- RsquareAdj(rda_int)$adj.r.squared

# compute the fractions of adjusted variation by subtraction
trt_R2 <- round(fra_d_R2 - fra_a_R2, digits = 2)
sn_R2  <- round(fra_d_R2 - fra_b_R2, digits = 2)

# sanity check
int_R2 <- round(fra_d_R2 - fra_c_R2, digits = 2)

# summary ----------------------------------------------------------------------
summary(rda_trt_sn)

# test of significance ---------------------------------------------------------

# test of individual RDA axes
anova(rda_trt_sn,by="axis", step = 1000) 
