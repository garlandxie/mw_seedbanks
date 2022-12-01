# libraries --------------------------------------------------------------------
library(janitor) # for cleaning column names
library(here)    # for creating relative file paths
library(dplyr)   # for manipulating data
library(tidyr)   # for cleaning data 

# import data ------------------------------------------------------------------

## seed bank -------------------------------------------------------------------
sb_spr <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

sb_fall <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

## plant status ----------------------------------------------------------------

plant_status <- read.csv(
  here("data", "intermediate_data", "plant_status.csv"),
  row.names = 1
)

# clean data -------------------------------------------------------------------

# merge both fall and spring datasets
sb <- rbind(sb_fall, sb_spr)

## proportions in the seed bank ------------------------------------------------

# calculate total number of seedling emergents)
sb_total <- sb %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(total_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, total_abund)

# calculate total number of seedling emergents that have an exotic status
sb_exotic <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "E") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(spontan_exotic_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, spontan_exotic_abund)

# calculate total number of native seedling emergents that are not
# in the TRCA native seed mix (so, spontaneous species here)
sb_spontan <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "N") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(spontan_nativ_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, spontan_nativ_abund)

# calculate total number of seedling emergents that are 
# in the TRCA native seed mix
sb_seed_mix <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "SM") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(sm_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, sm_abund)

# calculate total number of seedling emergents that are invasive
sb_invasive <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "I") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(spontan_inv_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, spontan_inv_abund)

# create a multi-column key (for record-keeping purposes)
multi_key_id <- c("season", "site_code", "treatment", "plot")

# join all of the datasets to calculate proportions
props <- sb_total %>%
  
  full_join(sb_exotic,    by = multi_key_id) %>%
  full_join(sb_spontan,   by = multi_key_id) %>%
  full_join(sb_seed_mix,  by = multi_key_id) %>%
  full_join(sb_invasive,  by = multi_key_id) %>%
  
  mutate(
    spotan_exotic_abund    = tidyr::replace_na(spontan_exotic_abund, 0),
    spotan_nativ_abund     = tidyr::replace_na(spontan_nativ_abund, 0),
    sm_abund               = tidyr::replace_na(sm_abund, 0),
    spontan_invasive_abund = tidyr::replace_na(spontan_inv_abund, 0)
  ) %>%
  
  mutate(
    
    props_spotan   = spotan_nativ_abund/total_abund, 
    props_sm       = sm_abund/total_abund,
    props_exotic   = spontan_exotic_abund/total_abund,
    props_invasive = spontan_inv_abund/total_abund 
    
  ) %>%
  
  mutate(
    treatment = factor(
      treatment, 
      levels = c("RES", "TIL", "MOW")
    )
  )
