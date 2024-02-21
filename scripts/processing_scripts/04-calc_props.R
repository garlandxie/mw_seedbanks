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
    spontan_exotic_abund    = tidyr::replace_na(spontan_exotic_abund, 0),
    spontan_nativ_abund     = tidyr::replace_na(spontan_nativ_abund, 0),
    sm_abund                = tidyr::replace_na(sm_abund, 0),
    spontan_inv_abund       = tidyr::replace_na(spontan_inv_abund, 0)
  ) %>%
  
  mutate(
    
    props_spontan_native    = spontan_nativ_abund/total_abund, 
    props_sm                = sm_abund/total_abund,
    props_spontan_exotic    = spontan_exotic_abund/total_abund,
    props_spontan_invasive  = spontan_inv_abund/total_abund 
  
  )
## assign missing values as zeros ----------------------------------------------

# some final clean-up prior to regression modelling and data visualization
props_tidy <- props %>%
  mutate(
    props_spontan_native   = tidyr::replace_na(props_spontan_native, 0),
    props_sm               = tidyr::replace_na(props_sm, 0),
    props_spontan_exotic   = tidyr::replace_na(props_spontan_exotic, 0),
    props_spontan_invasive = tidyr::replace_na(props_spontan_invasive , 0) 
  ) %>%
  
  mutate(
    props_spontan_native   = round(props_spontan_native, digits = 2),
    props_sm               = round(props_sm, digits = 2),
    props_spontan_exotic   = round(props_spontan_exotic, digits = 2),
    props_spontan_invasive = round(props_spontan_invasive, digits = 2)
  )

# summary statistics -----------------------------------------------------------

sb_exotic <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "E")

## til ----------------------------------------------------------------------

sb_exotic_til <- sb_exotic %>%
  filter(treatment == "TIL" & exotic_status == "Yes") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

sb_spontan_til <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "TIL" & status == "N") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

sb_seed_mix_til <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "TIL" & status == "SM") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

sb_invasive_til <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "TIL" & status == "I") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

## mow -------------------------------------------------------------------------

sb_exotic_mow <- sb_exotic %>%
  filter(treatment == "MOW" & exotic_status == "Yes") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

sb_spontan_mow <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "MOW" & status == "N") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

sb_seed_mix_mow <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "MOW" & status == "SM") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

sb_invasive_mow <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "MOW" & status == "I") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

## seed mix --------------------------------------------------------------------

sb_seed_mix_res <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "RES" & status == "SM") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

sb_exotic_res <- sb_exotic %>%
  filter(treatment == "RES" & exotic_status == "Yes") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

sb_invasive_res <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "RES" & status == "I") %>%
  pull(spp_code) %>%
  unique() %>%
  length()

sb_spontan_res <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(treatment == "RES" & status == "N") %>%
  pull(spp_code) %>%
  unique() %>% 
  length()

# summary statistics -----------------------------------------------------------

multi_key_id <- c("season", "site_code", "treatment", "plot")

sm_seedlings <- sb_seed_mix %>%
  full_join(sb_total, by = multi_key_id) %>%
  mutate(sm_abund = replace(sm_abund, is.na(sm_abund), 0)) %>%
  group_by(treatment) %>%
  summarize(
    sm_abund = sum(sm_abund, na.rm = FALSE),
    total_abund = sum(total_abund, na.rm = FALSE)
    )

sm_species <- sb %>%
  left_join(plant_status, by = c("spp_code" = "code")) %>%
  filter(status == "SM") %>%
  group_by(treatment) %>%
  summarize(sr = length(unique(spp_code))) %>%
  ungroup()

# save to disk -----------------------------------------------------------------

write.csv(
  x = props_tidy, 
  file = here("data", "intermediate_data", "props.csv")
)  
