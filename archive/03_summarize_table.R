# libraries --------------------------------------------------------------------
library(here)   # for creating relative file-paths
library(dplyr)  # for manipulating data 
library(tidyr)  # for making long to wide tables

# import -----------------------------------------------------------------------

## seed banks ------------------------------------------------------------------
sb_spr <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

sb_fall <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

## taxonomy --------------------------------------------------------------------
sb_taxon <- read.csv(
  here("data", "input_data", 
       "seed_bank_taxonomy.csv")
)

## plant status ----------------------------------------------------------------

plant_status  <- read.csv(
  here("data", "intermediate_data", 
       "plant_status.csv")
)

# clean data -------------------------------------------------------------------

sb <- rbind(sb_spr, sb_fall)

abund_by_site <- sb %>%
  full_join(sb_taxon, by = c("spp_code" = "Code")) %>%
  filter(!(is.na(plot))) %>%
  select(season, site_code, plot, spp_code, Genus, Species, total_abund) %>%
  group_by(season, site_code, spp_code, Genus, Species) %>%
  summarize(abund = sum(total_abund)) %>%
  ungroup()

abund_by_site_long <- abund_by_site %>%
  group_by(season, site_code, spp_code, Genus, Species) %>%
  tidyr::pivot_wider(names_from = site_code, values_from = abund) %>%
  ungroup() %>%
  mutate(across(M4_1:T2_3, ~replace_na(.x ,0)))

final_df <- abund_by_site_long %>%
  inner_join(plant_status, by = c("spp_code" = "code")) %>%
  select(season, spp_code, Genus, Species, M4_1:T2_3, status) %>%
  arrange(desc(season), spp_code) 
  
# save to disk -----------------------------------------------------------------

write.csv(
  x = final_df, 
  file = here("output", "data_appendix_output", "table_s2_summary_table.csv"), 
  row.names = FALSE
)


