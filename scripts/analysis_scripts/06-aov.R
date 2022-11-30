# libraries ----
library(here)        # for creating relative file-paths
library(ggplot2)     # for visualizing data 
library(dplyr)       # for manipulating data 
library(ggrepel)     # for creating easy-to-read text labels in ggplot2
library(tidyr)       # for creating wide or long data frames
library(vegan)       # for performing non-metric dimensional scaling
library(tibble)      # for assigning row names for tibbles
library(lme4)        # for running GLMM models
library(performance) # for running regression model diagnostics
library(emmeans)     # for calculating pairwise comparisons 
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

# clean data -------------------------------------------------------------------

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
  mutate(status = "Exotic") %>%
  dplyr::select(species = scientific_name, status)
  
## invasive status -------------------------------------------------------------

# obtain list of fifty invasive species in Toronto, Canada
invasive_plants_tidy <- invasive_plants %>%
  mutate(status = "Invasive") %>%
  dplyr::select(status, species = Species)

## seed mix status -------------------------------------------------------------

# obtain list of native plants across all curated seed mixes
# curated by the Toronto Regional Conservation Authority
# https://themeadoway.ca/resources/
seed_mix_tidy <- seed_mix %>%
  
  # clean columns to an R-friendly format
  janitor::clean_names() %>%
  
  # get binomial latin name (species and genus)
  mutate(species = paste(genus, species, sep = " ")) %>%
  
  # assign "seed mix" status to all species
  mutate(status = "Seed Mix") %>%
  
  # select relevant columns
  dplyr::select(species, status)

# taxonomy ----
sb_taxon_tidy <- sb_taxon %>%
  janitor::clean_names() %>%
  mutate(binom_latin = paste(genus, species, sep = " ")) %>%
  left_join(plants_to, by = c("binom_latin" = "SCIENTIFIC_NAME")) %>%
  dplyr::select(code, binom_latin, status = EXOTIC_NATIVE) %>%
  
  # manually assign exotic/native status 
  mutate(status = case_when(
    binom_latin == "Coreopsis tripteris" ~ "N",
    binom_latin == "Panicum milliaceum" ~ "E",
    binom_latin == "Picris echioides" ~ "E",
    binom_latin == "Verbenum urticifolia" ~ "N",
    binom_latin == "Veronica agrestis" ~ "E",
    binom_latin == "Euphorbia serpyllifolia" ~ "N", 
    binom_latin == "Acalypha rhomboidea" ~ "N",
    binom_latin == "Conyza canadensis" ~ "N",
    binom_latin == "Lactuca sativa" ~ "E",
    binom_latin == "Mentha arvensis" ~ "E",
    binom_latin == "Pycnanthemum muticum" ~ "N",
    binom_latin == "Sinapis alba" ~ "E",
    binom_latin == "Solidago rigida" ~ "N",
    TRUE ~ status)
  ) %>%
  
  # clean seed mix status
  left_join(seed_mix_S4, by = "binom_latin") %>%
  mutate(seed_mix_1_and_2 = replace_na(seed_mix_1_and_2, "No")) %>%
  select(code, binom_latin, status, seed_mix_1_and_2)

## proportions in the seed bank ----

# merge both fall and spring datasets
sb <- rbind(sb_fall, sb_spr)

# calculate the community abundance (i.e., total number of seedling germinants) 
# surveyed in each plot
sb_total <- sb %>%
  janitor::clean_names() %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(total_abund= sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, total_abund)

# calculate the total number of seedling germinants that have an exotic status
sb_exotic <- sb %>%
  janitor::clean_names() %>%
  left_join(sb_taxon_tidy, by = c("spp_code" = "code")) %>%
  filter(status == "E") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(exotic_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, exotic_abund)

# calculate the total number of seedling germinants that are not
# in the TRCA native seed mix (so, spontaneous species here)
sb_spotan <- sb %>%
  janitor::clean_names() %>%
  left_join(sb_taxon_tidy, by = c("spp_code" = "code")) %>%
  filter(seed_mix_1_and_2 == "No") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(spotan_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, spotan_abund)

# calculate the total number of seedling germinants that are 
# in the TRCA native seed mix
sb_seed_mix <- sb %>%
  janitor::clean_names() %>%
  left_join(sb_taxon_tidy, by = c("spp_code" = "code")) %>%
  filter(seed_mix_1_and_2 == "Yes") %>%
  group_by(season, site_code, treatment, plot) %>%
  summarize(sm_abund = sum(total_abund)) %>%
  ungroup() %>%
  dplyr::select(season, treatment, site_code, plot, sm_abund)

# create a multi-column key (for record-keeping purposes)
multi_key_id <- c("season", "site_code", "treatment", "plot")

# join all of the datasets to calculate proportions
props <- sb_total %>%
  
  full_join(sb_exotic,   by = multi_key_id) %>%
  full_join(sb_spotan,   by = multi_key_id) %>%
  full_join(sb_seed_mix, by = multi_key_id) %>%
  
  mutate(
    exotic_abund = tidyr::replace_na(exotic_abund, 0),
    spotan_abund = tidyr::replace_na(spotan_abund, 0),
    sm_abund     = tidyr::replace_na(sm_abund, 0)
    
  ) %>%
    
  mutate(
    
    props_spotan = spotan_abund/total_abund, 
    props_sm     = sm_abund/total_abund,
    props_exotic = exotic_abund/total_abund
    
  ) %>%
    
  mutate(
    treatment = factor(
      treatment, 
      levels = c("RES", "TIL", "MOW")
    )
  )

# plots ----

props_data_viz  <- props %>%
  mutate(
    treatment = as.character(treatment),
    treatment = case_when(
      treatment == "MOW" ~ "Maintenance-Mow", 
      treatment == "TIL" ~ "Tilling", 
      treatment == "RES" ~ "Undisturbed", 
      TRUE ~ treatment 
    ) 
  ) %>%
  mutate(season = factor(season, levels = c("Spring", "Fall")))

## proportion of exotics ----

(prop_exotic_plot <- props_data_viz %>%
  ggplot(aes(x = treatment, y = props_exotic, fill = treatment)) +
  geom_boxplot() + 
  geom_point(alpha = 0.5) + 
  labs(x = NULL, y = "Proportion of Exotic Species") + 
  theme_bw() + 
  facet_wrap(~season) + 
  theme(legend.position = "none")
)

## proportion of spontaneous species ---- 
(prop_spotan_plot <- props_data_viz %>%
   ggplot(aes(x = treatment, y = props_spotan, fill = treatment)) +
   geom_boxplot() + 
   geom_point(alpha = 0.5) + 
   labs(x = NULL, y = "Proportion of Spotaneous Species") + 
   facet_wrap(~season) + 
   theme_bw() + 
   theme(legend.position = "none")
)

## proportion of seed mix species ----
(prop_sm_plot <- props_data_viz %>%
   ggplot(aes(x = treatment, y = props_sm, fill = treatment)) +
   geom_boxplot() + 
   geom_point(alpha = 0.5) + 
   facet_wrap(~season) + 
   ylim(0, 1) + 
   labs(x = NULL, y = "Proportion of Native Species from the Seed Mix") + 
   theme_bw() + 
   theme(legend.position = "none")
)
  
# regression models ----

## proportion of exotics ----
prop_exotic_lm <- props %>%
  mutate(props_exotic = car::logit(props_exotic)) %>%
  lmer(
  props_exotic ~ treatment + (1|site_code),
  data = .
  )

summary(prop_exotic_lm)
performance::check_model(prop_exotic_lm)
performance::r2(prop_exotic_lm)

## proportion of spontaneous species ----
props_spotan_lm <-  props %>%
  mutate(logit_props_spotan = car::logit(props_spotan)) %>%
  lmer(
  logit_props_spotan ~ treatment + (1|site_code),
  data = .
  )

summary(props_spotan_glmer)
performance::check_model(props_spotan_glmer)
performance::r2(props_spotan_lm)

## proportion of seed mix species ----
props_sm_lm <- props %>%
  mutate(logit_props_sm = car::logit(props_sm)) %>%
  lmer(
  props_sm ~ treatment + (1|site_code),
  data = .
  )

summary(props_sm_glmer)
performance::check_model(props_sm_glmer)
performance::r2(props_sm_lm)

# pairwise comparisons ----
exotic_emm_trt <- emmeans(
  prop_exotic_glmer, 
  "treatment", 
  lmer.df = "satterthwaite"
)

pairs_exotic_emm <- pairs(exotic_emm_trt)

write.csv(
  x = pairs_exotic_emm, 
  file = here("output", "results", "pairs_exotic_emm.csv")
)

## save to disk ----

ggsave(
  here("output", "results", "prop_exotic.svg"), 
  plot = prop_exotic_plot,
  device = "svg", 
  units = "in",
  height = 5, 
  width = 6
)

ggsave(
  here("output", "results", "prop_spotan.svg"), 
  plot = prop_spotan_plot,
  device = "svg", 
  units = "in",
  height = 5, 
  width = 6 
)

ggsave(
  here("output", "results", "prop_sm.svg"), 
  plot = prop_sm_plot,
  device = "svg", 
  units = "in",
  height = 5, 
  width = 6 
)
