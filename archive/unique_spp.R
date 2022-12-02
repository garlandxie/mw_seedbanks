# libraries ----
library(here)
library(dplyr)

# import ----

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

sb_taxon <- read.csv(
  here("data", "input_data", 
       "seed_bank_taxonomy.csv")
)

# clean data ----

sb <- rbind(spr_sb, fall_sb)

# find unique species ----

# obtain a vector of species that emerged across the management regimes
sb_til <- sb %>%
  dplyr::filter(treatment == "TIL") %>%
  pull(spp_code) %>%
  unique()

sb_mow <- sb %>%
  dplyr::filter(treatment == "MOW") %>%
  pull(spp_code) %>%
  unique()

sb_res <- sb %>%
  dplyr::filter(treatment == "RES") %>%
  pull(spp_code) %>%
  unique()

# find unique species in tilled ----

# what species are unique to tilled plots that is not in mowed plots
til_mow <- setdiff(sb_til, sb_mow)

# what species are unique to tilled plots that is not in undisturbed plots
unique_til <- setdiff(til_mow, sb_res)

# find unique species in mowed ----

# what species are unique to tilled plots that is not in mowed plots
mow_til <- setdiff(sb_mow, sb_til)

# what species are unique to tilled plots that is not in undisturbed plots
unique_mow <- setdiff(mow_til, sb_res)

# find unique species in undisturbed ----

# what species are unique to tilled plots that is not in mowed plots
res_til <- setdiff(sb_res, sb_til)

# what species are unique to tilled plots that is not in undisturbed plots
unique_res <- setdiff(res_til, sb_mow)

