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
library(emmeans)     
library(ggsignif)    # 

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

# import dataset using R DRYAD API  
dryad_doi <- "10.5061/dryad.1ns1rn8sg"
dryad_link <- rdryad::dryad_download(dryad_doi)
plants_to <- read.csv(unlist(dryad_link))

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

# clean data ----

## seed mix ----

seed_mix_S4 <- seed_mix %>%
  janitor::clean_names() %>%
  
  # get binomial latin name (species and genus)
  mutate(binom_latin = paste(genus, species, sep = " ")) %>%
  
  # get species that were both in Meadoway seed mix 1 and 2 
  filter(seed_mix_1 == "Yes" | seed_mix_2 == "Yes") %>%
  mutate(seed_mix_1_and_2 = case_when(
    seed_mix_1 == "Yes" | seed_mix_2 == "Yes" ~ "Yes"
  )
  ) %>%
  dplyr::select(binom_latin, seed_mix_1, seed_mix_2, seed_mix_1_and_2)

## taxonomy ----
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

# regression models ----

props <- props %>%
  mutate(logit_props_exotic = car::logit(props_sm)) 

prop_exotic_lm <- lmer(
  logit_props_exotic ~ treatment + season + (1|site_code), 
  data = props
)

# pairwise comparisons ----

# calculate estimated marginal means
exotic_emm_trt <- emmeans(
  prop_exotic_lm, 
  "treatment", 
  lmer.df = "satterthwaite"
)

pairs_exotic_trt <- as.data.frame(pairs(exotic_emm_trt))

# obtain p-values for comparison between
# undisturbed and tilled
pairs_til_res <- pairs_exotic_trt %>%
  dplyr::filter(contrast == "RES - TIL") %>%
  pull(p.value) %>%
  signif(digits = 2)

# obtain p-value for comparison between 
# maintenance-mowing and undisturbed
pairs_mow_res <- pairs_exotic_trt %>%
  filter(contrast == "RES - MOW") %>%
  pull(p.value) %>%
  signif(digits = 2)

# obtain p-value for comparison between
# maintenance-mowing and tilling
pairs_mow_til <- pairs_exotic_trt %>%
  filter(contrast == "TIL - MOW") %>%
  pull(p.value) %>%
  signif(digits = 2)

# visualize the data ----
(prop_exotic_plot <- props_data_viz %>%
    ggplot(aes(x = treatment, y = props_exotic, fill = season)) +
    geom_boxplot() + 
    ylim(0, 3) + 
   
    geom_signif(
      y_position = 1.5, 
      xmin = 2, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_mow_res)),
      alpha = 0.5
    ) +
    geom_signif(
      y_position = 1.3, 
      xmin = 1, 
      xmax = 2, 
      annotation = paste("p = ", as.character(pairs_til_res)),
      alpha = 0.5
    ) + 
    geom_signif(
      y_position = 1.1, 
      xmin = 1, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_mow_til)),
      alpha = 0.5
    ) + 
   
    labs(
      x = "Management Regime", 
      y = "Proportion of Exotic Species"
    ) + 
    scale_fill_discrete(name = "Season") + 
    scale_y_continuous(breaks = seq(0, 1.0, 0.25)) + 
    theme_bw() + 
    theme(
      axis.title.x = element_text(
        margin = margin(t = 10, r = 0, b = 0, l = 0)
      ),
      axis.title.y = element_text(
        margin = margin(t = 0, r = 10, b = 0, l = 0)
      ),
      text = element_text(size = 15)  
    )
)

# save to disk -----

ggsave(
  filename = here("output", "results", "figure-S2.svg"), 
  device = "svg", 
  units = "in",
  height = 5, 
  width = 7
)

