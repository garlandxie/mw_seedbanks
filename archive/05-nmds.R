# libraries ----
library(here)    # for creating relative file-paths
library(ggplot2) # for visualizing data 
library(dplyr)   # for manipulating data 
library(ggrepel) # for creating easy-to-read text labels in ggplot2
library(tidyr)   # for creating wide or long data frames
library(vegan)   # for performing non-metric dimensional scaling
library(tibble)  # for assigning row names for tibbles

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
     "meadoway_seed_bank_taxonomy.csv")
)

## seed mix ----

seed_mix <- read.csv(
  here("data", "input_data", 
       "meadoway_seed_mix.csv")
)


# clean data ----

## seed mix ----

seed_mix_S4 <- seed_mix %>%
  janitor::clean_names() %>%
  
  # get binomial latin name (species and genus)
  mutate(binom_latin = paste(genus, species, sep = " ")) %>%
  
  # get species that were both in Meadoway seed mix 1 and 2 
  filter(seed_mix_1 == "Yes" | seed_mix_2 == "Yes") %>%
  mutate(
    seed_mix_1_and_2 = case_when(
      seed_mix_1 == "Yes" | seed_mix_2 == "Yes" ~ "Yes"
    )
  ) %>%
  dplyr::select(binom_latin, seed_mix_1, seed_mix_2, seed_mix_1_and_2)

## taxonomy ----
sb_taxon_tidy <- sb_taxon %>%
  janitor::clean_names() %>%
  mutate(binom_latin = paste(genus, species, sep = " ")) %>%
  left_join(plants_to, by = c("binom_latin" = "SCIENTIFIC_NAME")) %>%
  left_join(seed_mix_S4, by = "binom_latin") %>%
  dplyr::select(code, binom_latin, status = EXOTIC_NATIVE, seed_mix_1_and_2) %>%
  
  # manually assign exotic/native status 
  mutate(status = case_when(
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
  
  # manually assign native species that are in the seed mix
  mutate(status = case_when(
    seed_mix_1_and_2 == "Yes" ~ "SM", 
    TRUE ~ status)
  ) %>%
  
  mutate(status = case_when(
    is.na(binom_latin) | is.na(status) ~ "U", 
    TRUE ~ status)
  )

## seed bank ----

sb_tidy <- sb_fall %>%
  rbind(sb_spr) %>%
  janitor::clean_names() %>%
  group_by(section, site, treatment, plot, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup()

sb_comm_matrix <- sb_tidy %>%
  group_by(section, site, treatment, plot) %>%
  pivot_wider(
    names_from = spp_code, values_from = abundance) %>%
  ungroup() %>%
  mutate(across(everything(.), replace_na, 0)) 

sb_nmds <- sb_comm_matrix %>%
  dplyr::select(-section, -treatment, -site, -plot) %>%
  metaMDS(distance = "bray", trymax = 50)

saveRDS(
  object = sb_nmds, 
  file = here("data", "intermediate_data", "sb_nmds.rds")
)


stressplot(sb_nmds)


# plotting ----

## nmds: seed bank ----

sb_plot_scores <- 
  sb_nmds %>%
  scores() %>%
  as.data.frame() %>%
  cbind(sb_comm_matrix %>% dplyr::select(section, treatment, site, plot)) %>%
  dplyr::select(section, treatment, site, plot, NMDS1, NMDS2) %>%
  filter(!(site == "BRIM" & plot == 13))

sb_spp_scores <- 
  sb_nmds %>%
  scores("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "code") %>%
  left_join(sb_taxon_tidy, by = "code") %>%
  filter(code != "ECCR")

sb_nmds_plot <- ggplot() + 
    geom_point(
      data = sb_spp_scores,
      aes(x = NMDS1, y = NMDS2, col = status)
    ) +  
    geom_point(
      data = sb_plot_scores, 
      aes(x = NMDS1, y = NMDS2,
          shape = treatment), 
      alpha = 0.2
    ) + 
    
    scale_color_discrete(
      name = "Status", 
      labels = c(
        "Exotic (Spontaneous)", 
        "Native (Spotaneous)", 
        "Native (Seed Mix)", 
        "Unidentified")
    ) + 
    
    scale_shape_manual(
      name = "Management Regime", 
      labels = c("Maintenance-Mowed", "Undisturbed", "Disturbed"), 
      values = c(0, 2, 5)
    ) + 
    
    theme_bw()

# save to disk -----

ggsave(
  filename = here("output", "results", "sb_nmds.svg"),
  device = "svg", 
  units = "in", 
  height = 5, 
  width = 7
)


