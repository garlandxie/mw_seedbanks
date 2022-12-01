# libraries ----
library(here)        # for creating relative file-paths
library(dplyr)       # for manipulating data 
library(readxl)      # for reading excel files 
library(stringr)     # for manipulating string characters

# import ----

## seed bank ----
sb_spr <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

sb_fall <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

## exotic/native status ----

# import dataset using R Dryad API 
# Cadotte. 2021. Ecological Solutions and Evidence
# https://doi.org/10.1002/2688-8319.12036
dryad_cadotte_doi <- "10.5061/dryad.1ns1rn8sg"
dryad_cadotte_link <- rdryad::dryad_download(dryad_cadotte_doi)
plants_to <- read.csv(unlist(dryad_cadotte_link))

## invasive status ----

# import dataset using R Dryad API
# Potgieter et al. 2022. Journal of Applied Ecology
#  https://doi.org/10.1111/1365-2664.14103
dryad_potgieter_doi <- "10.5061/dryad.h9w0vt4k3"
dryad_potgieter_link <- rdryad::dryad_download(dryad_potgieter_doi)
dryad_potgieter_excel <- unlist(dryad_potgieter_link)[2]

invasive_plants <- 
  readxl::read_excel(
    dryad_potgieter_excel, 
    sheet = "Combined Species Ranking"
    )

## taxonomy ----

sb_taxon <- read.csv(
  here("data", "input_data", 
       "seed_bank_taxonomy.csv")
)

## seed mix ----

seed_mix <- read.csv(
  here("data", "input_data", 
       "seed_mix.csv")
)

# clean data: obtain plant status ----------------------------------------------

## native_status ---------------------------------------------------------------

# obtain list of native plants in Toronto, Canada
# can include plants in the TRCA curated seed mix as well
native_plants_tidy <- plants_to %>%
  janitor::clean_names() %>%
  filter(exotic_native == "N") %>%
  mutate(native_status = "Yes") %>%
  dplyr::select(species = scientific_name, native_status)

## exotic status ---------------------------------------------------------------

# obtain list of exotic plants in Toronto, Canada
# can include invasive plants as well 
exotic_plants_tidy <- plants_to %>%
  janitor::clean_names() %>%
  filter(exotic_native == "E") %>%
  mutate(exotic_status = "Yes") %>%
  dplyr::select(species = scientific_name, exotic_status)
  
## invasive status -------------------------------------------------------------

# obtain list of fifty invasive species in Toronto, Canada
invasive_plants_tidy <- invasive_plants %>%
  mutate(invasive_status = "Yes") %>%
  dplyr::select(species = Species, invasive_status)

## seed mix status -------------------------------------------------------------

# obtain list of native plants across all curated seed mixes
# curated by the Toronto Regional Conservation Authority
# https://themeadoway.ca/resources/
seed_mix_tidy <- seed_mix %>%
    janitor::clean_names() %>%
    mutate(species = paste(genus, species, sep = " ")) %>%
    mutate(seed_mix_status = "Yes") %>%
    dplyr::select(species, seed_mix_status)

# clean data: taxonomy ---------------------------------------------------------

sb_taxon_tidy <- sb_taxon %>%
  janitor::clean_names() %>%
  mutate(binom_latin = paste(genus, species, sep = " ")) %>%
  dplyr::select(code, binom_latin, authority)

# clean data: join taxonomy with all plant statuses ----------------------------

# do a left join of all plant status to the taxonomy df
# plant status: native, exotic, invasive, and seed mix
# each column will have a plant status 
# for transparency, in case there's a mistake that needs to be corrected

all_status <- sb_taxon_tidy %>%
  left_join(native_plants_tidy,   by = c("binom_latin" = "species")) %>%
  left_join(exotic_plants_tidy,   by = c("binom_latin" = "species")) %>%
  left_join(invasive_plants_tidy, by = c("binom_latin" = "species")) %>%
  left_join(seed_mix_tidy,        by = c("binom_latin" = "species"))




#  left_join(plants_to, by = c("binom_latin" = "SCIENTIFIC_NAME")) %>%
#  dplyr::select(code, binom_latin, status = EXOTIC_NATIVE) %>%
  
  # manually assign exotic/native status 
#  mutate(status = case_when(
#    binom_latin == "Coreopsis tripteris" ~ "N",
#    binom_latin == "Panicum milliaceum" ~ "E",
#    binom_latin == "Picris echioides" ~ "E",
#    binom_latin == "Verbenum urticifolia" ~ "N",
#    binom_latin == "Veronica agrestis" ~ "E",
#    binom_latin == "Euphorbia serpyllifolia" ~ "N", 
#    binom_latin == "Acalypha rhomboidea" ~ "N",
#    binom_latin == "Conyza canadensis" ~ "N",
#    binom_latin == "Lactuca sativa" ~ "E",
#    binom_latin == "Mentha arvensis" ~ "E",
#    binom_latin == "Pycnanthemum muticum" ~ "N",
#    binom_latin == "Sinapis alba" ~ "E",
#    binom_latin == "Solidago rigida" ~ "N",
#    TRUE ~ status)
#  ) %>%
  
#  # clean seed mix status
#  left_join(seed_mix_S4, by = "binom_latin") %>%
#  mutate(seed_mix_1_and_2 = replace_na(seed_mix_1_and_2, "No")) %>%
#  select(code, binom_latin, status, seed_mix_1_and_2)


